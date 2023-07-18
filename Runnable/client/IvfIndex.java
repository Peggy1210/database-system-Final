package org.vanilladb.core.storage.index.ivf;

import java.sql.DatabaseMetaData;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.vanilladb.core.sql.Schema;
import org.vanilladb.core.sql.Type;
import org.vanilladb.core.sql.IntegerConstant;
import org.vanilladb.core.server.VanillaDb;
import org.vanilladb.core.sql.BigIntConstant;
import org.vanilladb.core.sql.Constant;
import org.vanilladb.core.sql.VectorConstant;
import org.vanilladb.core.sql.distfn.EuclideanFn;
import org.vanilladb.core.storage.buffer.Buffer;
import org.vanilladb.core.storage.file.BlockId;
import org.vanilladb.core.storage.index.Index;
import org.vanilladb.core.storage.index.SearchKey;
import org.vanilladb.core.storage.index.SearchKeyType;
import org.vanilladb.core.storage.index.SearchRange;
import org.vanilladb.core.storage.metadata.TableInfo;
import org.vanilladb.core.storage.metadata.index.IndexInfo;
import org.vanilladb.core.storage.record.RecordFile;
import org.vanilladb.core.storage.record.RecordId;
import org.vanilladb.core.storage.record.RecordPage;
import org.vanilladb.core.storage.tx.Transaction;
import org.vanilladb.core.util.CoreProperties;

public class IvfIndex extends Index {
    private VectorConstant[] centroids; // TODO: Store in record file
    private Long[] adjMap;
    
    /* ================================ K-means ============================ */
    /* clustering - return the int list consists of group index of each data */
    public List<Integer> getNewCluster(List<DataEntry> rawData, VectorConstant[] centroids, int size) {
        List<Integer> newCluster = new ArrayList<Integer>();
        for (int i = 0; i < size; i++) {
            double minDist = Double.MAX_VALUE;
            int minCluster = -1;
            for (int j = 0 ; j < NUM_CLUSTERS ; j++) {
                double dist = getDistance(rawData.get(i), centroids[j]);
                if (dist < minDist) {
                    minDist = dist;
                    minCluster = j;
                }
            }
            newCluster.add(minCluster);
        }
        return newCluster;
    }

    /* calculate new centroids */
    public VectorConstant[] getNewCentroid(List<Integer> group_list, List<DataEntry> rawData, VectorConstant[] centroids, int size){
        VectorConstant [] newCentroids = new VectorConstant[NUM_CLUSTERS];
        for (int i = 0 ; i < NUM_CLUSTERS ; i++) {
            VectorConstant temp = VectorConstant.zeros(NUM_DIMENSIONS);
            int cnt = 0;
            for (int j = 0 ; j < size ; j++) {
                if (group_list.get(j) == i) {
                    SearchKey key = rawData.get(j).getSearchKey();
                    VectorConstant vec = (VectorConstant) key.get(KEY_ID);
                    temp = (VectorConstant) temp.add(vec);
                    cnt++;
                }
            }
            newCentroids[i] = (VectorConstant) temp.div(cnt);
        }
        return newCentroids;
    }

    /* calculate dist of two vector */
    /*****  used by `cluster` *******/
    private double getDistance(DataEntry entry, VectorConstant vec) {
        SearchKey key = entry.getSearchKey();
        return getDistance((VectorConstant) key.get(KEY_ID), vec);
    }
    
    public double getDistance(VectorConstant v1, VectorConstant v2) {
        ii.fieldNames();
        EuclideanFn distFun = new EuclideanFn("dist");
        distFun.setQueryVector(v1);
        return distFun.distance(v2);
        // double sum = 0;
        // for (int i = 0; i < v1.dimension(); i++) {
        //     double diff = v1.get(i) - v2.get(i);
        //     sum += diff * diff;
        // }
        // return Math.sqrt(sum);
    }

    /* check the stop condition of k-means */
    public boolean compareCentroids(VectorConstant[] newCentroids, VectorConstant[] oldCentroids, double threshold) {
        boolean isSame = false;
        double vdiff = 0;
        for (int i = 0; i < NUM_CLUSTERS; i++) {
            vdiff += getDistance(newCentroids[i], oldCentroids[i]);
        }
        if (vdiff / NUM_CLUSTERS < threshold) isSame = true;
        return isSame;
    }
    /* ================================ K-means ============================ */
    
    /*
     * Field name of the schema for centroid files.
     */
    static final String SCH_VECTOR = "vector", SCH_NUM_DATA = "num_data";
    /*
     * Field names of the schema for data files.
     */ 
    static final String SCH_KEY = "key", SCH_RID_BLOCK = "block", SCH_RID_ID = "id";
    /*
     * Field names of the schema for update log files.
     */ 
    static final String SCH_NUM_UPDATE = "num_update";

    private static final double KMEANS_THRESHOLD;
    private static final int NUM_CLUSTERS, NUM_UPDATES;
    private static final int NUM_DIMENSIONS;
    private static final int CHOOSE_N; // Enable choosing N nearest cluster
    private static final int KEY_ID = 0;
    private static final int TOP_K;
    private RecordFile centroidFile;
    private RecordFile dataFile;
    private RecordFile updateLogFile;
    private int update_cnt;

    static {
        NUM_CLUSTERS = CoreProperties.getLoader().getPropertyAsInteger(
                IvfIndex.class.getName() + ".NUM_CLUSTERS", 20);
        NUM_UPDATES = CoreProperties.getLoader().getPropertyAsInteger(
            IvfIndex.class.getName() + ".NUM_UPDATES", 20);
        KMEANS_THRESHOLD =  CoreProperties.getLoader().getPropertyAsDouble(
            IvfIndex.class.getName() + ".KMEANS_THRESHOLD", 1e2);
        NUM_DIMENSIONS = CoreProperties.getLoader().getPropertyAsInteger(
                IvfIndex.class.getName() + ".NUM_DIMENSIONS", 48);
        CHOOSE_N = CoreProperties.getLoader().getPropertyAsInteger(
            IvfIndex.class.getName() + ".CHOOSE_N", 10);
        TOP_K = CoreProperties.getLoader().getPropertyAsInteger(
            IvfIndex.class.getName() + ".TOP_K", 20);
    }
    
    private int searchCnt = 0;
    private int[] candidateCentroid;
    private int[] numData;
    private boolean isBeforeFirsted;
    private boolean clustered;
    private boolean isUpdate;
    private List<DataEntry> sortedData;

    public IvfIndex(IndexInfo ii, SearchKeyType keyType, Transaction tx) {
        super(ii, keyType, tx);
        centroids = new VectorConstant[NUM_CLUSTERS];
        numData = new int[NUM_CLUSTERS];
        candidateCentroid = new int[CHOOSE_N];
        getCentroidFile();
        isUpdate = false;
        sortedData = new ArrayList<DataEntry>();
        // System.out.println("Constructor");
    }

    public static long searchCost(SearchKeyType keyType, long totRecs, long matchRecs) {
    	// Calculate search cost
        int rpb = Buffer.BUFFER_SIZE / RecordPage.slotSize(clusterSchema(keyType));
        return (totRecs / rpb) / NUM_CLUSTERS;
        // return 1;
    }

    private static String keyFieldName(int index) {
        return SCH_KEY + index;
    }

    private static Schema clusterSchema(SearchKeyType keyType) {
        Schema sch = new Schema();
        for (int i = 0; i < keyType.length(); i++)
            sch.addField(keyFieldName(i), keyType.get(i));
        sch.addField(SCH_RID_BLOCK, Type.BIGINT);
        sch.addField(SCH_RID_ID, Type.INTEGER);
        return sch;
    }

    private static Schema centroidSchema() {
        Schema sch = new Schema();
        sch.addField(SCH_VECTOR, Type.VECTOR(NUM_DIMENSIONS));
        sch.addField(SCH_NUM_DATA, Type.INTEGER);
        return sch;
    }

    private static Schema updateLogSchema() {
        Schema sch = new Schema();
        sch.addField(SCH_NUM_UPDATE, Type.INTEGER);
        return sch;
    }

    @Override
    public void beforeFirst(SearchRange searchRange) {
        // throw new UnsupportedOperationException("Unimplemented method 'beforeFirst'");
        if (!searchRange.isSingleValue())
            throw new UnsupportedOperationException();
        
        SearchKey searchKey = searchRange.asSearchKey();
        beforeFirst(searchKey);
    }

    public void beforeFirst(SearchKey searchKey) {
        close();

        if (clustered) {
            // If the data has been classified, find the N nearest centroids
            findCandidateCentroid(searchKey);

            if (isUpdate) {
                // Doing updates, read the data file corresponding to the nearest cluster
                String tblName = ii.indexName() + "-cluster" + candidateCentroid[0];
                TableInfo ti = new TableInfo(tblName, clusterSchema(keyType));
                dataFile = ti.open(tx, false);

                if (dataFile.fileSize() == 0) RecordFile.formatFileHeader(ti.fileName(), tx);
                dataFile.beforeFirst();
                // System.out.println("go to " + candidateCentroid[0]);
            } else {
                printDataFile();
                List<DataEntry> data = new ArrayList<DataEntry>();
                for (int i = 0; i < CHOOSE_N; i++) {
                    data.addAll(getDataFile(candidateCentroid[i]));
                }
                int dataSize = data.size();
                // System.out.println("total " + dataSize + " found");

                // Calculate the distance of each data with respect to the target
                Map<Double, DataEntry> distMap = new HashMap<Double, DataEntry>();
                List<Double> dist = new ArrayList<>();
                VectorConstant v = (VectorConstant) searchKey.get(KEY_ID);
                for (int i = 0; i < dataSize; i++) {
                    double d = getDistance(v, (VectorConstant) data.get(i).getSearchKey().get(KEY_ID));
                    distMap.put(d, data.get(i));
                    dist.add(d);
                }

                // Sort and put into new data list
                dist.sort(Comparator.naturalOrder());

                // Put data into `sortedData`
                sortedData.clear();
                for (int i = 0; i < dataSize; i++) sortedData.add(distMap.get(dist.get(i)));
                // System.out.println("total " + sortedData.size() + " found");
            }
        } else {
            String tblName = ii.indexName() + "-cluster0";
            TableInfo ti = new TableInfo(tblName, clusterSchema(keyType));
            dataFile = ti.open(tx, false);

            if (dataFile.fileSize() == 0) RecordFile.formatFileHeader(ti.fileName(), tx);
            dataFile.beforeFirst();
        }
        isBeforeFirsted = true;
        searchCnt = 0;
    }

    private void findCandidateCentroid(SearchKey searchKey) {
        // TODO: Find nearest clusters
        VectorConstant v = (VectorConstant) searchKey.get(KEY_ID);

        // Calculate the nearest cluster with respect to the given vector
        // XXX: use adjacency map of centroids
        Map<Double, Integer> distMap = new HashMap<Double, Integer>();
        List<Double> dist = new ArrayList<>();
        for (int i = 0; i < NUM_CLUSTERS; i++) {
            double d = getDistance(v, centroids[i]);
            distMap.put(d, i);
            dist.add(d);
        }

        // sort
        dist.sort(Comparator.naturalOrder());

        // find the order of candidate centroid
        for (int i = 0; i < CHOOSE_N; i++) {
            candidateCentroid[i] = distMap.get(dist.get(i));
        }
    }

    @Override
    public boolean next() {
        if (!isBeforeFirsted)
            throw new IllegalStateException("Must call beforeFirst() before calling next()");

        if (isUpdate || !clustered) {
            // If it is updating data, call next() of the opened data file
            // Also, if the data is not clustered, we will use `dataFile` to stored data
            return dataFile.next();
        } else {
            // If it is searching, return the sorted data list
            searchCnt++;
            // System.out.println(sortedData.get(searchCnt));
            if (searchCnt == sortedData.size() -1) return false;
            else return true;
        }
    }

    /*
     * Close all the opened files
     */
    @Override
    public void close() {
        centroidFile.close();
        if (dataFile != null) dataFile.close();
    }

    @Override
    public RecordId getDataRecordId() {
        if (isUpdate || !clustered) {
            // If it is updating data, return the RecordId of the opened data file
            // If the data has not been clustered, use `dataFile`
            long blkNum = (Long) dataFile.getVal(SCH_RID_BLOCK).asJavaVal();
            int id = (Integer) dataFile.getVal(SCH_RID_ID).asJavaVal();
            return new RecordId(new BlockId(dataFileName, blkNum), id);
        } else {
            // If it is searching, return the current RecordId of the sorted data
            long blkNum = sortedData.get(searchCnt).getBlkNum();
            int id = sortedData.get(searchCnt).getId();
            RecordId rid = new RecordId(new BlockId(dataFileName, blkNum), id);
            // System.out.println(rid);
            return rid;
            // return new RecordId(new BlockId(dataFileName, blkNum), id);
        }
    }

    @Override
    public void insert(SearchKey key, RecordId dataRecordId, boolean doLogicalLogging) {
        isUpdate = true;
        beforeFirst(key);
        
        // log logical operation start
        if (doLogicalLogging)
            tx.recoveryMgr().logLogicalStart();

        dataFile.insert();
        for (int i = 0; i < keyType.length(); i++)
            dataFile.setVal(keyFieldName(i), key.get(i));
        dataFile.setVal(SCH_RID_BLOCK, new BigIntConstant(dataRecordId.block().number()));
        dataFile.setVal(SCH_RID_ID, new IntegerConstant(dataRecordId.id()));

        checkUpdates();

        // log logical operation end
        if (doLogicalLogging)
            tx.recoveryMgr().logIndexInsertionEnd(ii.indexName(), key,
                    dataRecordId.block().number(), dataRecordId.id());
        isUpdate = false;
    }

    @Override
    public void delete(SearchKey key, RecordId dataRecordId, boolean doLogicalLogging) {
        isUpdate = true;
        beforeFirst(key);

        // log logical operation start
        if (doLogicalLogging)
            tx.recoveryMgr().logLogicalStart();

        while (next())
            if (getDataRecordId().equals(dataRecordId)) {
                dataFile.delete();
                return;
            }

        checkUpdates();

        // log logical operation end
        if (doLogicalLogging)
            tx.recoveryMgr().logIndexDeletionEnd(ii.indexName(), key,
                    dataRecordId.block().number(), dataRecordId.id());

        isUpdate = false;
    }

    private long fileSize(String fileName) {
		tx.concurrencyMgr().readFile(fileName);
		return VanillaDb.fileMgr().size(fileName);
	}

    public void checkUpdates() {
        // Update num_update in the centroid file
        getUpdateLogFile();
        update_cnt++;
        if (update_cnt >= NUM_UPDATES) {
            runKMeans();
            resetUpdateLogFile();
            clustered = true;
        } else {
            setUpdateLogFile(update_cnt);
        }
        updateLogFile.close();
    }

    @Override
    public void preLoadToMemory() {
        String tblName = ii.indexName() + "-centroid.tbl";
        long size = fileSize(tblName);
        BlockId blk;
        for (int j = 0; j < size; j++) {
            blk = new BlockId(tblName, j);
            tx.bufferMgr().pin(blk);
        }
    }

    private void getUpdateLogFile() {
        String tblName = ii.indexName() + "-update-log";
        TableInfo ti = new TableInfo(tblName, updateLogSchema());
        updateLogFile = ti.open(tx, false);
        if (centroidFile.fileSize() == 0) RecordFile.formatFileHeader(ti.fileName(), tx);
        
        updateLogFile.beforeFirst();
        if(updateLogFile.next())
            update_cnt = (Integer) updateLogFile.getVal(SCH_NUM_UPDATE).asJavaVal();
        else 
            update_cnt = 0;
    }

    private void resetUpdateLogFile() {
        setUpdateLogFile(0);
    }

    private void setUpdateLogFile(int updates) {
        String tblName = ii.indexName() + "-update-log";
        VanillaDb.fileMgr().delete(tblName + ".tbl");        
        TableInfo ti = new TableInfo(tblName, updateLogSchema());
        updateLogFile = ti.open(tx, false);
        if (updateLogFile.fileSize() == 0) RecordFile.formatFileHeader(ti.fileName(), tx);
        updateLogFile.beforeFirst();
        updateLogFile.insert();
        updateLogFile.setVal(SCH_NUM_UPDATE, new IntegerConstant(updates));
    }

    private boolean set = false;
    public void setCentroidFile(Map<Integer, List<DataEntry>> dataList) {
        // Construct a new centroid file
        String tblName = ii.indexName() + "-centroid";
        VanillaDb.fileMgr().delete(tblName + ".tbl");
        TableInfo ti = new TableInfo(tblName, centroidSchema());
        centroidFile = ti.open(tx, false);
        if (centroidFile.fileSize() == 0) RecordFile.formatFileHeader(ti.fileName(), tx);
        centroidFile.beforeFirst();
        // int cnt = 0;
        for (int c = 0; c < NUM_CLUSTERS; c++) {
            centroidFile.insert();
            centroidFile.setVal(SCH_VECTOR, centroids[c]);
            // cnt += dataList.get(c).size();
            centroidFile.setVal(SCH_NUM_DATA, new IntegerConstant(dataList.get(c).size()));
        }
        // System.out.println(cnt);
        // if (cnt == 10000) set = true;
    }

    public void getCentroidFile() {
        // Load the centroid information
        String tblName = ii.indexName() + "-centroid";
        TableInfo ti = new TableInfo(tblName, centroidSchema());
        centroidFile = ti.open(tx, false);
        if (centroidFile.fileSize() == 0) RecordFile.formatFileHeader(ti.fileName(), tx);

        centroidFile.beforeFirst();
        if (centroidFile.fileSize() == 1) {
            // There is only file header in the centroid file
            clustered = false;
        } else {
            int id = 0;
            centroidFile.beforeFirst();
            while(centroidFile.next()) {
                centroids[id] = (VectorConstant) centroidFile.getVal(SCH_VECTOR);
                numData[id] = (Integer) centroidFile.getVal(SCH_NUM_DATA).asJavaVal();
                id++;
            }
            clustered = true;
        }
    }

    public void setDataFile(List<DataEntry> data, int clusterId) {
        String tblName = ii.indexName() + "-cluster" + clusterId;
        VanillaDb.fileMgr().delete(tblName + ".tbl");
        TableInfo ti = new TableInfo(tblName, clusterSchema(keyType));
        RecordFile tempDataFile = ti.open(tx, false);
        if (tempDataFile.fileSize() == 0) RecordFile.formatFileHeader(ti.fileName(), tx);
        
        // if (set && clusterId == 157) System.out.println("populating cluster " + clusterId);
        tempDataFile.beforeFirst();
        for (DataEntry d: data) {
            SearchKey key = d.getSearchKey();
            Long blk = d.getBlkNum();
            int id = d.getId();
            tempDataFile.insert();
            for (int i = 0; i < keyType.length(); i++)
                tempDataFile.setVal(keyFieldName(i), key.get(i));
            tempDataFile.setVal(SCH_RID_BLOCK, new BigIntConstant(blk));
            tempDataFile.setVal(SCH_RID_ID, new IntegerConstant(id));
            // if (set && clusterId == 157) System.out.print("["+blk+","+id+"]");
        }
        // if (set && clusterId == 157) System.out.print("\n");

        // int cnt = 0;
        // if (set && clusterId == 157) {
        //     tempDataFile.beforeFirst();
        //     while(tempDataFile.next()) {
        //         // cnt++;
        //         Long blk = (Long) tempDataFile.getVal(SCH_RID_BLOCK).asJavaVal();
        //         int id = (Integer) tempDataFile.getVal(SCH_RID_ID).asJavaVal();
        //         System.out.print("["+blk+","+id+"]");
        //     }
        //     System.out.print("\n");
        // }
        
        // Close the file
        tempDataFile.close();
    }

    public List<DataEntry> getDataFile(int clusterId) {
        RecordFile tempDataFiles;
        // Open data file
        String tblName = ii.indexName() + "-cluster" + clusterId;
        TableInfo ti = new TableInfo(tblName, clusterSchema(keyType));
        tempDataFiles = ti.open(tx, false);
        if (tempDataFiles.fileSize() == 0) RecordFile.formatFileHeader(ti.fileName(), tx);
        tempDataFiles.beforeFirst();

        // Extract the information of the data file
        List<DataEntry> data = new ArrayList<DataEntry>();
        while(tempDataFiles.next()) {
            // Extract the data information
            long blkNum = (Long) tempDataFiles.getVal(SCH_RID_BLOCK).asJavaVal();
            int id = (Integer) tempDataFiles.getVal(SCH_RID_ID).asJavaVal();
            Constant[] val = new Constant[keyType.length()];
            for (int j = 0; j < keyType.length(); j++)
                val[j] = (Constant) tempDataFiles.getVal(keyFieldName(j));
            DataEntry e = new DataEntry(blkNum, id, val);
            data.add(e);

            // if (set && clusterId == 157) System.out.print("["+blkNum+","+id+"]");
        }
        // if (set && clusterId == 157) System.out.print("\n");
        // Close the file
        tempDataFiles.close();
        return data;
    }

    // Map<DataEntry, Integer> map;

    public void printDataFile() {
        // int cnt = 0;
        // for (int i = 0; i < NUM_CLUSTERS; i++) {
            // System.out.println("cluster " + i);
            List<DataEntry> list = getDataFile(157);
            for (int j = 0; j < list.size(); j++){
                System.out.print("["+list.get(j).getBlkNum()+","+list.get(j).getId()+"]");
            }
            System.out.println("\ntotal: " + list.size());
            // cnt += list.size();
        // }
        // System.out.println("All data: " + cnt);
    }

    public void runKMeans() {
        // if there are more than NUM_UPDATE updates, run KMeans to update centroid
        // load data
        Map<Integer, List<DataEntry>> dataList = new HashMap<Integer, List<DataEntry>>();
        for (int i = 0; i < NUM_CLUSTERS; i++) {
            List<DataEntry> data = getDataFile(i);
            dataList.put(i, data);
        }
        
        // prepare data
        int totalDataCount = 0;
        for (int i = 0 ; i < NUM_CLUSTERS; i++) totalDataCount += dataList.get(i).size();
        List<Integer> group_list = new ArrayList<Integer>(); 
        List<DataEntry> rawData = new ArrayList<DataEntry>();
        VectorConstant[] newCentroids = new VectorConstant[NUM_CLUSTERS];
        
        if (dataList.get(1).size() == 0) { /* First time simulation */
            rawData = dataList.get(0);
            for (int i = 0; i < NUM_CLUSTERS; i++) {
                DataEntry e = dataList.get(0).get(i);
                centroids[i] = (VectorConstant) e.getSearchKey().get(KEY_ID);
            }
        } else {
            /* Update */
            for (int c = 0; c < NUM_CLUSTERS; c++)
                rawData.addAll(dataList.get(c));
        }

        // run KMeans
        while(true){
            group_list = getNewCluster(rawData, centroids, totalDataCount);
            newCentroids = getNewCentroid(group_list, rawData, centroids, totalDataCount);
            if (compareCentroids(newCentroids, centroids, KMEANS_THRESHOLD)) break;
            else centroids = Arrays.copyOf(newCentroids, newCentroids.length);
        }

        Map<Integer, List<DataEntry>> newDataList = new HashMap<Integer, List<DataEntry>>();
        for (int i = 0; i < NUM_CLUSTERS; i++) newDataList.put(i, new ArrayList<DataEntry>());
        for (int i = 0; i < totalDataCount; i++) 
            newDataList.get(group_list.get(i)).add(rawData.get(i));

        // store new clusters
        setCentroidFile(newDataList);
        for (int c = 0; c < NUM_CLUSTERS; c++) {
            List<DataEntry> data = newDataList.get(c);
            setDataFile(data, c);
        }
    }
}
