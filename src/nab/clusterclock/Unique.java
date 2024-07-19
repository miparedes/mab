package nab.clusterclock;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.core.BEASTObject;
import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Tree;


@Description("calculates sum of a valuable")
public class Unique extends CalculationNode implements Function, Loggable {
    final public Input<List<Function>> functionInput = new Input<>("arg", "argument to be summed", new ArrayList<>(), Validate.REQUIRED);


    enum Mode {integer_mode, double_mode}

    Mode mode;

    boolean needsRecompute = true;
    double unique = 0;
    double storedUnique = 0;

    @Override
    public void initAndValidate() {
        List<Function> valuable = functionInput.get();
        mode = Mode.integer_mode;
        for (Function v : valuable) {
	        if (!(v instanceof IntegerParameter || v instanceof BooleanParameter)) {
	            mode = Mode.double_mode;
	        }
        }
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return unique;
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
    	unique = 1;
    	List<Double> values = new ArrayList<>();
        for (Function v : functionInput.get()) {
            for (int i = 0; i < v.getDimension(); i++) {
            	values.add(v.getArrayValue(i));
            }
        }
        Collections.sort(values);
        for (int i = 1; i < values.size(); i++) {
        	if (Math.abs(values.get(i-1)-values.get(i))>0.5)
                unique++;
        }
        needsRecompute = false;
    }

    @Override
    public double getArrayValue(int dim) {
        if (dim == 0) {
            return getArrayValue();
        }
        return Double.NaN;
    }

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
        storedUnique = unique;
        super.store();
    }

    @Override
    public void restore() {
    	unique = storedUnique;
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }

    /**
     * Loggable interface implementation follows
     */
    @Override
    public void init(PrintStream out) {
        out.print("unique(" + ((BEASTObject) functionInput.get().get(0)).getID() + ")\t");
    }

    @Override
    public void log(long sampleNr, PrintStream out) {
        double sum = 1;
    	List<Double> values = new ArrayList<>();
        for (Function v : functionInput.get()) {
            for (int i = 0; i < v.getDimension(); i++) {
            	values.add(v.getArrayValue(i));
            }
        }
        Collections.sort(values);
        for (int i = 1; i < values.size(); i++) {
        	if (Math.abs(values.get(i-1)-values.get(i))>0.5)
        		sum++;
        }

        out.print((int) sum + "\t");
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

} // class Sum
