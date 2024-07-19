package nab.util;


import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.core.*;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;

import org.apache.commons.math.MathException;


@Description("Produces prior (log) probability of value x." +
        "If x is multidimensional, the components of x are assumed to be independent, " +
        "so the sum of log probabilities of all elements of x is returned as the prior.")
public class LargerThanPrior extends Distribution {
    final public Input<Function> m_x = new Input<>("x", "point at which the density is calculated", Validate.REQUIRED);

    /**
     * shadows distInput *
     */
    protected ParametricDistribution dist;

    @Override
    public void initAndValidate() {
        calculateLogP();
    }

    @Override
    public double calculateLogP() {
    	logP=0;
        Function x = m_x.get();
        if (x instanceof RealParameter || x instanceof IntegerParameter) {
            // test that parameter is inside its bounds
            double value = x.getArrayValue(0);
            double value2 = x.getArrayValue(1);
            if (value>value2) {
                logP = Double.NEGATIVE_INFINITY;
                return Double.NEGATIVE_INFINITY;
            }
            
        }
        return logP;
    }

    /**
     * return name of the parameter this prior is applied to *
     */
    public String getParameterName() {
        if (m_x.get() instanceof BEASTObject) {
            return ((BEASTObject) m_x.get()).getID();
        }
        return m_x.get() + "";
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


}
