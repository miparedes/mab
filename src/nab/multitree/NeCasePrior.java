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
public class NeCasePrior extends CalculationNode implements Function, Loggable {
    final public Input<Function> functionInput = new Input<>("arg", "argument for which the differences for entries is calculated", Validate.REQUIRED);
    final public Input<Function> casesInput = new Input<>("logCases", "log of the cases", Validate.REQUIRED);
    final public Input<Function> overallNeScalerInput = new Input<>("overallNeScaler", "argument for which the differences for entries is calculated", Validate.REQUIRED);

    enum Mode {integer_mode, double_mode}

    Mode mode;

    boolean needsRecompute = true;
    double[] errorTerm;
    double[] storedErrorTerm;

    @Override
    public void initAndValidate() {
    	errorTerm = new double[functionInput.get().getDimension()];
    	storedErrorTerm = new double[functionInput.get().getDimension()];
    }

    @Override
    public int getDimension() {
        return errorTerm.length;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return errorTerm[0];
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
    	
		for (int i = 0; i < functionInput.get().getDimension(); i++) {

			errorTerm[i] = functionInput.get().getArrayValue(i) - casesInput.get().getArrayValue(i)- overallNeScalerInput.get().getArrayValue(0);
//			System.out.println(functionInput.get().getDimension());
//			System.out.println(casesInput.get().getDimension());
//			System.out.println(overallNeScalerInput.get().getDimension());
//			System.out.println(errorTerm.length);			
		}
        needsRecompute = false;
    }

    @Override
    public double getArrayValue(int dim) {
        if (needsRecompute) {
            compute();
        }
        return errorTerm[dim];
    }

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
    	System.arraycopy(errorTerm, 0, storedErrorTerm, 0, errorTerm.length);
        super.store();
    }

    @Override
    public void restore() {
    	double [] tmp = storedErrorTerm;
    	storedErrorTerm = errorTerm;
    	errorTerm = tmp;
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }

	@Override
	public void init(PrintStream out) {
		for (int i = 1; i < errorTerm.length; i++) {
			out.print("errorTerm_"+i+"\t");
			
		}
		
	}
	
	@Override
	public void log(long iteration, PrintStream out) {
		for (int i = 1; i < errorTerm.length; i++) {
			out.print(errorTerm[i] + "\t");
		}
				
		
	}

	@Override
	public void close(PrintStream out) {
	}
} // class Sum