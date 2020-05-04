package nab.skygrid;


import java.util.Collections;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;


/**
 * @author Nicola F. Mueller
 */
@Description("Rate with time vector.")
public class TimeVaryingRates extends CalculationNode {
	
    final public Input<RealParameter> rateInput = new Input<>("rate",
            "Nes over time in log space", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts",
            "When to switch between elements of Ne", Input.Validate.REQUIRED);

    RealParameter rate;
    RealParameter rateShifts;
    
    double[] growth;
    double[] growth_stored;


    @Override
	public void initAndValidate() {
    	rate = rateInput.get();
    	rateShifts = rateShiftsInput.get();
    	rate.setDimension(rateShifts.getDimension());  
    	recalculateRate();
    	
//    	System.out.print("rate=[");
//    	double time = 0.0;
//    	while (time < 0.15) {
//    		System.out.print(getMeanRate(time, time+0.001) + ",");
//    		time+=0.0001;
//    	}
//    	System.out.print("];");
//    	System.exit(0);
    			
    }

	public double getRate(double t) {		
		
		int intervalnr = getIntervalNr(t);
		if (intervalnr>=rateShifts.getDimension()) {
			return Double.POSITIVE_INFINITY;
		}		
		
		double timediff = t;
		if (intervalnr>0)
			timediff -= rateShifts.getArrayValue(intervalnr-1);

//		return rate.getArrayValue(intervalnr)-growth[intervalnr]*timediff;
		return rate.getArrayValue(intervalnr);
	}
	
	public double getMeanRate(double start, double finish) {
    	if (start==finish)
    		return 0.0;
    	
    	// get the interval "start" is in
    	int first_int = getIntervalNr(start);
    	// get the interval "finish" is in
    	int last_int = getIntervalNr(finish);
    	
    	double weighted = 0.0;
    	double curr_time = start; 
    	double tot_time = 0.0;

    	for (int i = first_int; i <= last_int;i++) {
    		if (i>=rateShifts.getDimension()) {
    			return Double.POSITIVE_INFINITY;
    		}    		
    		double next_time = Math.min(getTime(i), finish);
    		
    		double timediff1 = curr_time;
    		double timediff2 = next_time;
    		if (i>0) {
    			timediff1 -= rateShifts.getArrayValue(i-1);
    			timediff2 -= rateShifts.getArrayValue(i-1);
    		}

//    		weighted +=  (timediff2-timediff1)*((rate.getArrayValue(i) - growth[i]*timediff2) + (rate.getArrayValue(i) - growth[i]*timediff1));

    		weighted +=  (next_time - curr_time)*rate.getArrayValue(i);
    		curr_time = next_time;
    		tot_time += (timediff2-timediff1);
    	}
    	return weighted/tot_time;
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
	
	// computes the Ne's at the break points
	private void recalculateRate() {
		growth = new double[rateShifts.getDimension()];
		double curr_time = 0.0;
		for (int i = 1; i < rate.getDimension(); i++) {
			growth[i-1] = (rate.getArrayValue(i-1)- rate.getArrayValue(i))/(rateShifts.getArrayValue(i-1)-curr_time);
			curr_time = rateShifts.getArrayValue(i-1);
		}
	}
	
	@Override
	public boolean requiresRecalculation() {
		recalculateRate();
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