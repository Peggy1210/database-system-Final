package org.vanilladb.bench.server.param.ann;

import org.vanilladb.core.sql.Schema;
import org.vanilladb.core.sql.storedprocedure.SpResultRecord;
import org.vanilladb.core.sql.storedprocedure.StoredProcedureHelper;

public class AnnTestbedLoaderParamHelper implements StoredProcedureHelper {

    private static final String TABLES_DDL[] = new String[1];
    private static final String INDEX_DDL[] = new String[1];

    private int numOfItems, numDimension;

    public String[] getTableSchemas() {
        return TABLES_DDL;
    }

    public int getNumberOfItems() {
        return numOfItems;
    }

    public int getVecDimension() {
        return numDimension;
    }

    // TODO:
    public String[] getIndexSchemas() {
        return INDEX_DDL;
    }

    @Override
    public void prepareParameters(Object... pars) {
        numOfItems = (Integer) pars[0];
        numDimension = (Integer) pars[1];
        TABLES_DDL[0] = "CREATE TABLE items (i_id INT, i_emb VECTOR(" + numDimension + "), i_name VARCHAR(24))";
        INDEX_DDL[0] = "CREATE INDEX idx_item ON items (i_emb) USING IVF";
    }

    @Override
    public Schema getResultSetSchema() {
        return new Schema();
    }

    @Override
    public SpResultRecord newResultSetRecord() {
        return new SpResultRecord();
    }

    @Override
    public boolean isReadOnly() {
        return false;
    }
}
