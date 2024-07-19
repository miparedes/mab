package nab.clusterclock;


import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.core.*;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;


@Description("Produces prior (log) probability of value x." +
        "If x is multidimensional, the components of x are assumed to be independent, " +
        "so the sum of log probabilities of all elements of x is returned as the prior.")
public class PrePrior extends Distribution {
    final public Input<IntegerParameter> m_x = new Input<>("x", "point at which the density is calculated", Validate.REQUIRED);
    final public Input<RealParameter> pdfInput = new Input<>("pdf", "input of pdf values for int parameters, start counting at 0.", Validate.REQUIRED);

    /**
     * shadows distInput *
     */

    @Override
    public void initAndValidate() {
        calculateLogP();
    }

    @Override
    public double calculateLogP() {
    	logP = 0.0;
    	for (int i = 0; i < m_x.get().getDimension(); i++) {
    		logP += Math.log(pdfInput.get().getArrayValue((int) m_x.get().getArrayValue(i)));
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
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        conditions.add(pdfInput.get().getID());
        return conditions;
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();

        String id = null;
        if (m_x.get() != null && m_x.get() instanceof BEASTInterface) {
            arguments.add(((BEASTInterface)m_x.get()).getID());
        }

        return arguments;
    }

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
}
