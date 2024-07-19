package epochclock;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.PopulationFunction;

public class EpochClockRate extends BranchRateModel.Base {

	final public Input<PopulationFunction> timevaryingClockRatesInput = new Input<>("timevaryingClockRates",
			"the Clock rates through time.", Input.Validate.REQUIRED);

	@Override
	public void initAndValidate() {
	}
	

	@Override
	public double getRateForBranch(Node node) {			
		if (node.isRoot()) {
			return 0.0;		
		}	
		return timevaryingClockRatesInput.get().getIntegral(node.getHeight(), node.getParent().getHeight())/(node.getParent().getHeight()-node.getHeight());
	}


	@Override
	protected boolean requiresRecalculation() {
		if (((CalculationNode) timevaryingClockRatesInput.get()).isDirtyCalculation()) {
			// this is only called if any of its inputs is dirty, hence we need to recompute
			return true;
		}

		return false;
	}


}
