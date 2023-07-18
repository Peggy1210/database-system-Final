package org.vanilladb.core.sql.distfn;

import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

import org.vanilladb.core.sql.VectorConstant;
import org.vanilladb.core.util.CoreProperties;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;

public class EuclideanFn extends DistanceFn {

    private static final int USE_SIMD;

    static {
        USE_SIMD = CoreProperties.getLoader().getPropertyAsInteger(
                EuclideanFn.class.getName() + ".USE_SIMD", 1);
    }

    public EuclideanFn(String fld) {
        super(fld);
    }

    // DONE: SIMD version 2-norm
    @Override
    protected double calculateDistance(VectorConstant vec) {

        if (USE_SIMD == 1) {
            int[] v1 = query.asJavaVal();
            int[] v2 = vec.asJavaVal();
            
            VectorSpecies<Integer> species = IntVector.SPECIES_PREFERRED;
            double sumSqrDiff = 0;
            
            IntVector fv1, fv2, fv3;
            for (int i = 0 ; i < v1.length; i += species.length()) {
                VectorMask<Integer> m = species.indexInRange(i, vec.copy().length);
                
                fv1 = IntVector.fromArray(species, v1, i, m);
                fv2 = IntVector.fromArray(species, v2, i, m);
                fv3 = fv1.sub(fv2);

                // For some unknown reason, fv3.mul(fv3) is significantly faster than fv3.pow(2).
                sumSqrDiff += fv3.mul(fv3).reduceLanes(VectorOperators.ADD);
            }
            return Math.sqrt(sumSqrDiff);
        } else {
            double sum = 0;
            for (int i = 0; i < vec.dimension(); i++) {
                double diff = query.get(i) - vec.get(i);
                sum += diff * diff;
            }
            return Math.sqrt(sum);
        }   	
    }
}
