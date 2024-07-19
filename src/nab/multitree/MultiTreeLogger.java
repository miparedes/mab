package nab.multitree;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import beast.core.CalculationNode;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.IntervalType;
import beast.mascot.distribution.MascotNative2;
import beast.mascot.dynamics.Dynamics;
import beast.mascot.ode.*;
import cern.colt.Arrays;
import nab.skygrid.TimeVaryingRates;

/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a beast.tree using under the framework of Mueller (2017).")
@Citation("Nicola F. Müller, David A. Rasmussen, Tanja Stadler (2017)\n  The Structured Coalescent and its Approximations.\n  Mol Biol Evol 2017 msx186. doi: 10.1093/molbev/msx186")
public class MultiTreeLogger extends MultiTreeDistribution implements Loggable {
	
    final public Input<Integer> minClusterSizeInput = new Input<>("minClusterSize", 
    		"A population size model", 0);

	
	MultiTreeIntervals treeIntervals; 
    @Override
    public void initAndValidate() {
		treeIntervals = multiTreeIntervalsInput.get();
	}

    /**
     * Loggable interface implementation follows *
     */
    @Override
    public void init(final PrintStream out) {
        out.println("#NEXUS\n");
        out.println("Begin trees;\n");        
    }

    @Override
    public void log(final long sample, final PrintStream out) {
        calculateLogP();
        double maxHeight=-1.0;
        for (int i = 0; i < treeIntervals.treeInput.get().size();i++) {
        	maxHeight = Math.max(treeIntervals.treeInput.get().get(i).getRoot().getHeight()+
        			treeIntervals.offset[i] + treeIntervals.rootLengthInput.get().get(i).getValue(), 
        			maxHeight);
        }
        String tree_string = "rem";
        
        for (int i = 0; i < treeIntervals.treeInput.get().size();i++) {
        	if (treeIntervals.treeInput.get().get(i).getExternalNodes().size()>=minClusterSizeInput.get())
	        	tree_string = tree_string + "," + 
		        getNewick(treeIntervals.treeInput.get().get(i), i) + 
		        ":" +  (maxHeight-treeIntervals.treeInput.get().get(i).getRoot().getHeight()-treeIntervals.offset[i]-treeIntervals.rootLengthInput.get().get(i).getValue());
        }
        
        tree_string = tree_string.replace("rem,", "(");
        tree_string = tree_string + "):0.0";
        
        
//        System.out.println(Arrays.toString(nodeProbs.get(96)));
       
        
        out.print("tree STATE_" + sample + " = ");
        out.print(tree_string);
        out.print(";");	        		
    }

    private String getNewick(Tree tree, int i) {
        final StringBuilder rootString = new StringBuilder();

        rootString.append("(" + toNewick(tree.getRoot()));
        rootString.append(":" + treeIntervals.rootLengthInput.get().get(i).getValue());
        rootString.append(")[&obs=0]");
//        rootString.append(")" + getStateProbs(rootProbs.get(tree.getRoot().getNr()+nodeOffset)));
		return rootString.toString();
	}
    
    public String toNewick(Node n) {
        final StringBuilder buf = new StringBuilder();
        if (!n.isLeaf()) {
            buf.append("(");
            boolean isFirst = true;
            for (Node child : n.getChildren()) {
                if (isFirst)
                    isFirst = false;
                else
                    buf.append(",");
                
                buf.append(toNewick(child));
            }
            buf.append(")[&obs=1]");

            if (n.getID() != null)
                buf.append(n.getID());
        } else {
            buf.append(n.getID() + "[&obs=1]");
        }
        if (!n.isRoot())
        	buf.append(":").append(n.getLength());
        return buf.toString();
    }
    
	@Override
    public void close(final PrintStream out) {
        out.print("End;");
    }

}
