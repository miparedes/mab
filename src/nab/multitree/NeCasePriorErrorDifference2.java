package nab.multitree;


import java.io.PrintStream;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import beast.core.Function;
import beast.core.Loggable;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.CalculationNode;
import beast.core.Description;


@Description("calculates the differences between the entries of a vector")
public class NeCasePriorErrorDifference2 extends CalculationNode implements Function, Loggable {
    final public Input<Function> functionInput = new Input<>("arg", "argument for which the differences for entries is calculated", Validate.REQUIRED);
    final public Input<Function> casesInput = new Input<>("logCases", "log of the cases", Validate.REQUIRED);
    final public Input<Function> overallNeScalerInput = new Input<>("overallNeScaler", "argument for which the differences for entries is calculated", Validate.REQUIRED);
    final public Input<Function> rateShifts = new Input<>("rateShifts", "argument for which the differences for entries is calculated", Validate.REQUIRED);

    enum Mode {integer_mode, double_mode}

    Mode mode;

    boolean needsRecompute = true;
    double[] errorTerm;
    double[] storedErrorTerm;
    double[] errorTermDiff;
    double[] storedErrorTermDiff;
    

    @Override
    public void initAndValidate() {
    	errorTerm = new double[functionInput.get().getDimension()];
    	storedErrorTerm = new double[functionInput.get().getDimension()];
    	errorTermDiff = new double[functionInput.get().getDimension()-1];
    	storedErrorTermDiff = new double[functionInput.get().getDimension()-1];
    	
    }

    @Override
    public int getDimension() {
        return errorTermDiff.length;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return errorTermDiff[0];
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
    	
		for (int i = 0; i < functionInput.get().getDimension(); i++) {

			errorTerm[i] = functionInput.get().getArrayValue(i) - casesInput.get().getArrayValue(i)- overallNeScalerInput.get().getArrayValue(0);
			
		
		}

		double[] growthRates = new double[errorTerm.length-1];
        //loop over all time points
    	for (int j = 1; j < errorTerm.length; j++){
    		double errorTermDiff_1 = errorTerm[j] - errorTerm[j-1];
    		
    		double timediff = rateShifts.get().getArrayValue(j-1);
    		if (j>1)
    			timediff -= rateShifts.get().getArrayValue(j-2);

//    		double timediff = rateShifts.get().getArrayValue(j) - rateShifts.get().getArrayValue(j-1);
    		growthRates[j-1] = errorTermDiff_1/timediff;
    	}

//		System.out.println(functionInput.get().getDimension());
//		System.out.println(casesInput.get().getDimension());
//		System.out.println(overallNeScalerInput.get().getDimension());
//		System.out.println(errorTerm.length);	
//		System.out.println(growthRates.length);
//
//		System.out.println(errorTermDiff.length);
//		System.out.println(rateShifts.get().getDimension());
    	for (int j = 1; j < growthRates.length; j++){
    		 errorTermDiff[j-1] = growthRates[j]-growthRates[j-1];
    	}
    	
        // add contribution from first or last entry
    	errorTermDiff[errorTermDiff.length-1] = growthRates[growthRates.length-1];
        needsRecompute = false;
    }

    @Override
    public double getArrayValue(int dim) {
        if (needsRecompute) {
            compute();
        }
        return errorTermDiff[dim];
    }

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
    	System.arraycopy(errorTermDiff, 0, storedErrorTermDiff, 0, errorTermDiff.length);
        super.store();
    }

    @Override
    public void restore() {
    	double [] tmp = storedErrorTermDiff;
    	storedErrorTermDiff = errorTermDiff;
    	errorTermDiff = tmp;
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }

	@Override
	public void init(PrintStream out) {
		for (int i = 1; i < errorTermDiff.length; i++) {
			out.print("errorTermDiff_"+i+"\t");
			
		}
		
	}
	
	@Override
	public void log(long iteration, PrintStream out) {
		for (int i = 1; i < errorTermDiff.length; i++) {
			out.print(errorTermDiff[i] + "\t");
		}
				
		
	}

	@Override
	public void close(PrintStream out) {
	}
} // class Sum