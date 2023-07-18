package org.vanilladb.core.storage.index.ivf;

import org.vanilladb.core.sql.Constant;
import org.vanilladb.core.storage.index.SearchKey;

public class DataEntry {
    private SearchKey key;
    private Long blk;
    private int id;

    public DataEntry(Long blk, int id, SearchKey key) {
        this.blk = blk;
        this.id = id;
        this.key = key;
    }

    public DataEntry(Long blk, int id, Constant... val) {
        this.key = new SearchKey(val);
        this.blk = blk;
        this.id = id;
    }

    public SearchKey getSearchKey() {
        return key;
    }

    public Long getBlkNum() {
        return blk;
    }

    public int getId() {
        return id;
    }
}
