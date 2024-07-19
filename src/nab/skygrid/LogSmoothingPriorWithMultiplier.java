package nab.skygrid;

import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;

public class LogSmoothingPriorWithMultiplier extends Distribution {
	
    
    final public Input<ParametricDistribution> distInput = new Input<>("distr", 
    		"distribution used to calculate prior on the difference between intervals, e.g. normal, beta, gamma.", 
    		Validate.REQUIRED);
    
    final public Input<ParametricDistribution> initDistrInput = new Input<>("initialDistr", 
    		"distribution used to calculate prior on the difference between intervals, e.g. normal, beta, gamma.");
        
    final public Input<ParametricDistribution> finalDistrInput = new Input<>("finalDistr", 
    		"distribution used to calculate prior on the difference between intervals, e.g. normal, beta, gamma.", 
    		Input.Validate.XOR, initDistrInput);
    
    final public Input<RateMultiplier> rateMultiplierInput = new Input<>("rateMultiplier", 
    		"A population size model", Validate.REQUIRED);
    
    public Input<Skygrid> skygridInput = new Input<>(
    		"skygrid", "input of skygrid effective population sizes", Validate.REQUIRED);        


   
    
    private Skygrid skygrid;
    private RateMultiplier mult;
    
    protected ParametricDistribution dist;
    protected ParametricDistribution initDistr;
    protected ParametricDistribution finalDistr;
    
    @Override
    public void initAndValidate() {
    	skygrid = skygridInput.get();    	
        dist = distInput.get();
        mult = rateMultiplierInput.get();
        
        if (initDistrInput.get()!=null)
        	initDistr = initDistrInput.get();
        if (finalDistrInput.get()!=null)
        	finalDistr = finalDistrInput.get();
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
        double offset=0.0001;
        //loop over all time points
    	for (int j = 1; j < skygrid.rateShifts.getDimension(); j++){
    		double old_time = skygrid.rateShifts.getArrayValue(j-1)-offset;
			double new_time = skygrid.rateShifts.getArrayValue(j-1)-offset;
			
			double Ne_old = skygrid.getPopSize(old_time) + mult.getRate(old_time);
			double Ne_new = skygrid.getPopSize(new_time) + mult.getRate(new_time);
			
    		double diff = Ne_new-Ne_old;    		
    		logP += dist.logDensity(diff);
    	}
        
        // add contribution from first or last entry
        if (initDistrInput.get()!=null) {
			double Ne = skygrid.getPopSize(0.0) + mult.getRate(0.0);
    		logP += initDistr.logDensity(Ne);
        }
        
        return logP;
    }
}
