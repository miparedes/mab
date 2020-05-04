package nab.skygrid;

import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;


/**
 * @author Nicola F. Mueller
 */
@Description("Populaiton function with values at certain time points that are interpolated in between. Parameter has to be in log space")
public class Skygrowth extends PopulationFunction.Abstract {
	
    final public Input<RealParameter> NeInput = new Input<>("logNe",
            "Nes over time in log space", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts",
            "When to switch between elements of Ne", Input.Validate.REQUIRED);

    //
    // Public stuff
    //
    RealParameter Ne;
    RealParameter rateShifts;
    
    boolean NesKnown = false;
    double[] growth;
    double[] growth_stored;


    @Override
	public void initAndValidate() {
    	Ne = NeInput.get();    	    	
    	rateShifts = rateShiftsInput.get();
    	Ne.setDimension(rateShifts.getDimension()+1);
    	growth = new double[rateShifts.getDimension()];
    	recalculateNe();
    }


	@Override
	public List<String> getParameterIds() {
		throw new IllegalArgumentException("Not implemented");
	}


	@Override
	public double getPopSize(double t) {	
//		if (!NesKnown)
//			recalculateNe();

		int intervalnr = getIntervalNr(t);
		if (intervalnr>=rateShifts.getDimension()) {
			return Double.POSITIVE_INFINITY;
		}
		double timediff = t;
		if (intervalnr>0)
			timediff -= rateShifts.getArrayValue(intervalnr-1);
				
		return Math.exp(Ne.getArrayValue(intervalnr)-growth[intervalnr]*timediff);
	}

	
    @Override
	public double getIntegral(double start, double finish) {
//		if (!NesKnown)
//			recalculateNe();

    	if (start==finish)
    		return 0.0;
    	
    	// get the interval "start" is in
    	int first_int = getIntervalNr(start);
    	// get the interval "finish" is in
    	int last_int = getIntervalNr(finish);
    	
    	double weighted = 0.0;
    	double curr_time = start; 
    	
    	for (int i = first_int; i <= last_int;i++) {
    		if (i>=rateShifts.getDimension()) {
    			return Double.POSITIVE_INFINITY;
    		}
    		
    		double next_time = Math.min(getTime(i), finish);
    		double r = growth[i];  
    		
    		double timediff1 = curr_time;
    		double timediff2 = next_time;
    		if (i>0) {
    			timediff1 -= rateShifts.getArrayValue(i-1);
    			timediff2 -= rateShifts.getArrayValue(i-1);
    		}

    		if (r==0.0)
    			weighted +=  (next_time-curr_time)/Math.exp(Ne.getArrayValue(i));
    		else
    			weighted +=  (Math.exp(timediff2*r) - Math.exp(timediff1*r))/Math.exp(Ne.getArrayValue(i))/r;
    				
    		curr_time = next_time;
    	}
    	return weighted;
    }


	@Override
	public double getInverseIntensity(double x) {
		throw new IllegalArgumentException("Not implemented");
	}
	
	private int getIntervalNr(double t) {
		// check which interval t + offset is in
		for (int i = 0; i < rateShifts.getDimension(); i++)
			if (t<rateShifts.getArrayValue(i))
				return i;
		
		// after the last interval, just keep using the last element
		return rateShifts.getDimension();					
	}

	private double getTime(int i) {
		return rateShifts.getArrayValue(i);
	}

	@Override
	public double getIntensity(double t) {
		throw new IllegalArgumentException("Not implemented");
	}
	
	// computes the Ne's at the break points
	private void recalculateNe() {
		growth = new double[rateShifts.getDimension()];
		double curr_time = 0.0;
		for (int i = 1; i < Ne.getDimension(); i++) {
			growth[i-1] = (Ne.getArrayValue(i-1)- Ne.getArrayValue(i))/(rateShifts.getArrayValue(i-1)-curr_time);
			curr_time = rateShifts.getArrayValue(i-1);
		}
		NesKnown = true;
	}

	@Override
	public boolean requiresRecalculation() {
		recalculateNe();
		return super.requiresRecalculation();
	}
	
	@Override
	public void store() {
		growth_stored = new double[growth.length];
		System.arraycopy(growth, 0, growth_stored, 0, growth.length);
		super.store();
	}
	
	@Override
	public void restore() {
		System.arraycopy(growth_stored, 0, growth, 0, growth_stored.length);
		super.restore();
	}


	
}