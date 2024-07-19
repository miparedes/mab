package nab.skygrid;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;

public class LogSmoothingPriorIndependent extends Distribution {
	
    public Input<RealParameter> NeLogInput = new Input<>(
    		"NeLog", "input of effective population sizes");        
    
    final public Input<ParametricDistribution> distInput = new Input<>("distr", 
    		"distribution used to calculate prior on the difference between intervals, e.g. normal, beta, gamma.", 
    		Validate.REQUIRED);
    
    final public Input<ParametricDistribution> initDistrInput = new Input<>("initialDistr", 
    		"distribution used to calculate prior on the difference between intervals, e.g. normal, beta, gamma.",
    		Input.Validate.OPTIONAL);
        
    final public Input<ParametricDistribution> finalDistrInput = new Input<>("finalDistr", 
    		"distribution used to calculate prior on the difference between intervals, e.g. normal, beta, gamma.", 
    		Input.Validate.OPTIONAL);
    
    final public Input<List<ParametricDistribution>> jumpDistrInput = new Input<>("jumpDistr", 
    		"distribution used to calculate prior on the difference between intervals, e.g. normal, beta, gamma.", 
    		new ArrayList());
    
    final public Input<IntegerParameter> jumpsInput = new Input<>("jumps", 
    		"intervals when jumps are allowed.", 
    		Input.Validate.REQUIRED);


   
    
    private RealParameter NeLog;
    
    protected ParametricDistribution dist;
    protected ParametricDistribution initDistr;
    protected ParametricDistribution finalDistr;
    
    private List<Integer> jumpIntervals;
    
    @Override
    public void initAndValidate() {
    	NeLog = NeLogInput.get();    	
        dist = distInput.get();
        
        if (initDistrInput.get()!=null)
        	initDistr = initDistrInput.get();
        if (finalDistrInput.get()!=null)
        	finalDistr = finalDistrInput.get();

        jumpIntervals = new ArrayList<>();
        for (int i = 0; i < jumpsInput.get().getDimension(); i++) {
        	jumpIntervals.add((int) jumpsInput.get().getArrayValue(i));
        }
        
    }


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
	
    public double calculateLogP() {
    	logP = 0;
               
        
        //loop over all time points
    	for (int j = 1; j < NeLog.getDimension(); j++){
    		double diff = NeLog.getArrayValue(j) - NeLog.getArrayValue(j-1);    	
    		if (jumpIntervals.contains(j)) {
    			logP += jumpDistrInput.get().get(jumpIntervals.indexOf(j)).logDensity(diff);
    		}else {
    			logP += dist.logDensity(diff);
    		}
    	}
        
        // add contribution from first or last entry
        if (initDistrInput.get()!=null)
    		logP += initDistr.logDensity(NeLog.getArrayValue(0));
        if (finalDistrInput.get()!=null) {
    		logP += finalDistr.logDensity(NeLog.getArrayValue(NeLog.getDimension()-1));
        }       
        return logP;
    }
}
