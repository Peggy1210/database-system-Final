package org.vanilladb.core.sql;

/*
 * A ConstantRange for VectorConstant type. Currently VectorConstantRange does not support range.
 */
public class VectorConstantRange extends ConstantRange {

    private VectorConstant low, high;
    private boolean lowIncl, highIncl;

    public VectorConstantRange(VectorConstant low, boolean lowIncl, VectorConstant high, boolean highIncl) {
        // if (low.equals(high))
        //     throw new UnsupportedOperationException("VectorConstantRange does not support range");
        this.low = low;
        this.high = high;
        this.lowIncl = lowIncl;
        this.highIncl = highIncl;
    }

    @Override
    public boolean isValid() {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'isValid'");
    }

    @Override
    public boolean hasLowerBound() {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'hasLowerBound'");
    }

    @Override
    public boolean hasUpperBound() {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'hasUpperBound'");
    }

    @Override
    public Constant low() {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'low'");
    }

    @Override
    public Constant high() {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'high'");
    }

    @Override
    public boolean isLowInclusive() {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'isLowInclusive'");
    }

    @Override
    public boolean isHighInclusive() {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'isHighInclusive'");
    }

    @Override
    public double length() {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'length'");
    }

    @Override
    public ConstantRange applyLow(Constant c, boolean inclusive) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'applyLow'");
    }

    @Override
    public ConstantRange applyHigh(Constant c, boolean incl) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'applyHigh'");
    }

    @Override
    public ConstantRange applyConstant(Constant c) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'applyConstant'");
    }

    @Override
    public boolean isConstant() {
        // TODO Auto-generated method stub
        // throw new UnsupportedOperationException("Unimplemented method 'isConstant'");
        return true;
    }

    @Override
    public Constant asConstant() {
        // TODO Auto-generated method stub
        // throw new UnsupportedOperationException("Unimplemented method 'asConstant'");
        return low;
    }

    @Override
    public boolean contains(Constant c) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'contains'");
    }

    @Override
    public boolean lessThan(Constant c) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'lessThan'");
    }

    @Override
    public boolean largerThan(Constant c) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'largerThan'");
    }

    @Override
    public boolean isOverlapping(ConstantRange r) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'isOverlapping'");
    }

    @Override
    public boolean contains(ConstantRange r) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'contains'");
    }

    @Override
    public ConstantRange intersect(ConstantRange r) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'intersect'");
    }

    @Override
    public ConstantRange union(ConstantRange r) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'union'");
    }
    
}
