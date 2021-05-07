package bdmm.stateclock;

import java.io.PrintStream;
import java.util.Arrays;

import bdmmprime.mapping.TypeMappedTree;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;

public class BdmmStateClock extends BranchRateModel.Base implements Loggable {

	final public Input<RealParameter> relativeClockRatesInput = new Input<>("relativeClockRates",
			"the rate parameters associated with rates in particular states.", Input.Validate.REQUIRED);

	final public Input<FlatTypeMappedTree> typedTreeInput = new Input<>("typedTree", "bdmm typed mapped tre.",
			Input.Validate.REQUIRED);

	public int states;
	RealParameter relativeClockRates;
	FlatTypeMappedTree typedTree;
	
	boolean isMapped = false;

	@Override
	public void initAndValidate() {
		typedTree = typedTreeInput.get();

		relativeClockRates = relativeClockRatesInput.get();

		// ensure that the rate parameters has the same dimension as there are states in
		// the stateClockInput
		states = typedTree.parameterizationInput.get().getNTypes();
		if (relativeClockRates.getDimension() != states)
			relativeClockRates.setDimension(states);
	}

	@Override
	public double getRateForBranch(Node node) {			
		if (node.isRoot())
			return 0.0;
		
		if (!isMapped) {
			typedTree.doStochasticMapping();
			isMapped = true;
		}
		
		double t = 0;
		for (int i = 0; i < states; i++) {
			t += typedTree.getNodeTime(node.getNr(), i);
		}
		
//		System.out.println("node");
//
//		System.out.println(typedTree.getNodeTime(node.getNr(), 0) + " " + typedTree.getNodeTime(node.getNr(), 1));

		if (t-node.getParent().getHeight()+node.getHeight()>0.000001) {
			System.out.println("err");
			System.out.println(t + " " + (node.getParent().getHeight()-node.getHeight()));
			System.out.println(node.getHeight());
			System.out.println(typedTree.treeInput.get());
			System.exit(0);
		}
		
		
		// compute the mean rate
		double rate = 0.0;
		for (int i = 0; i < states; i++) {
			rate += typedTree.getNodeTime(node.getNr(), i) * relativeClockRates.getArrayValue(i) ;
		}
//		System.out.println(rate);
//		System.out.println(node.getParent().getHeight()+ " " + node.getHeight());
//		System.out.println(t);
//		System.out.println(rate * meanRateInput.get().getValue()/t);
		if (Double.isNaN(rate * meanRateInput.get().getValue()/t))
			System.out.println(t);
		return rate * meanRateInput.get().getValue()/t;
	}


	@Override
	protected boolean requiresRecalculation() {
		if (((CalculationNode) typedTree).isDirtyCalculation()) {
			// this is only called if any of its inputs is dirty, hence we need to recompute
			isMapped = false;
			return true;
		}
				
		if (meanRateInput.get().isDirty(0)) {
			isMapped = false;
			return true;
		}
		
		for (int i = 0; i < relativeClockRates.getDimension(); i++) {
			if (relativeClockRates.isDirty(i)) {
				isMapped = false;
				return true;
			}
		}
		
		if (super.requiresRecalculation()) {
			isMapped = false;
			return true;
		}
			

		return false;
	}

//	@Override
//	protected void store() {
//		super.store();
//	}
//
//	@Override
//	protected void restore() {
//		super.restore();
//	}

	@Override
	public void init(PrintStream out) {
		for (int i = 0; i < relativeClockRates.getDimension(); i++) {
			out.print(String.format("%s.%d\t", relativeClockRates.getID(), i));
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0; i < relativeClockRates.getDimension(); i++) {
			out.print(relativeClockRates.getArrayValue(i) + "\t");
		}
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub

	}

}
