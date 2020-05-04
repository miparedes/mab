package nab.multitree;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.tree.coalescent.TreeIntervals;

public class MultiTreeDistribution extends Distribution {
    final public Input<MultiTreeIntervals> multiTreeIntervalsInput = new Input<>("multiTreeIntervals", "Intervals for a phylogenetic beast tree", Validate.REQUIRED);

	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
    
    @Override
    protected boolean requiresRecalculation() {
    	MultiTreeIntervals ti = multiTreeIntervalsInput.get();
        if (ti != null) {
            //boolean d = ti.isDirtyCalculation();
            //assert d;
            assert ti.isDirtyCalculation();
            return true;
        }
    	
        return false;
    }


}
