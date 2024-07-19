package nab.skygrid;


import java.util.Collections;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;


/**
 * @author Nicola F. Mueller
 */
@Description("Rate with time vector.")
public class RateMultiplier extends CalculationNode {
	
    final public Input<RealParameter> casesInput = new Input<>("cases",
            "Number of cases from outside", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts",
            "When to switch between elements of Ne", Input.Validate.REQUIRED);
    final public Input<Double> minCasesInput = new Input<>("minCases",
            "minmal number of cases for log standardization", 1.0);


    RealParameter logStandardCases;
    RealParameter rateShifts;
    

    @Override
	public void initAndValidate() {
    	logStandardCases = casesInput.get();

    	rateShifts = rateShiftsInput.get();
    	if (logStandardCases.getDimension()!=rateShifts.getDimension())
    		throw new IllegalArgumentException("cases and rate shifts have different dimension");
    	// log everything
    	for (int i = 0; i < logStandardCases.getDimension(); i++)
    		logStandardCases.setValue(i, Math.log(Math.max(minCasesInput.get(),logStandardCases.getArrayValue(i))));
    	
    	System.out.println(logStandardCases);

    	// standardize everything
    	double mean=0.0, std=0.0;
    	for (int i = 0; i < logStandardCases.getDimension(); i++)
    		mean += logStandardCases.getArrayValue(i);
    	
    	mean/=logStandardCases.getDimension();

    	for (int i = 0; i < logStandardCases.getDimension(); i++)
    		std += Math.pow((logStandardCases.getArrayValue(i)-mean),2);
    	
    	std /= logStandardCases.getDimension()-1;
    	std = Math.sqrt(std);
    	
    	
    	for (int i = 0; i < logStandardCases.getDimension(); i++)
    		logStandardCases.setValue(i, (logStandardCases.getArrayValue(i)-mean)/std );

    }

	public double getRate(double t) {				
		int intervalnr = getIntervalNr(t);
		if (intervalnr>=rateShifts.getDimension()) {
			return Double.POSITIVE_INFINITY;
		}		
		return logStandardCases.getArrayValue(intervalnr);
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

    	for (int i = first_int; i <= last_int;i++) {
    		if (i>=rateShifts.getDimension()) {
    			return Double.POSITIVE_INFINITY;
    		}    		
    		double next_time = Math.min(getTime(i), finish);
    		weighted +=  (next_time - curr_time)*logStandardCases.getArrayValue(i);
    		curr_time = next_time;
    	}
    	return weighted/(finish-start);
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

}