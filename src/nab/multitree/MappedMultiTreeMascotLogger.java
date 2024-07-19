package nab.multitree;

import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;

public class MappedMultiTreeMascotLogger extends CalculationNode implements Loggable {
	
	final public Input<MappedMultitreeMascot> mappedMascotInput = new Input<>("mappedMultitreeMascot", 
			"tree for which to calculate the intervals", Input.Validate.REQUIRED);
	
    final public Input<Integer> minClusterSizeInput = new Input<>("minClusterSize", 
    		"A population size model", 0);

	MappedMultitreeMascot mmm;
	
	@Override
	public void initAndValidate() {
		mmm = mappedMascotInput.get();			
	}

	
	@Override
	public void init(PrintStream out) {
        out.println("#NEXUS\n");
        out.println("Begin trees;\n");        
	}	
	
	@Override
	public void log(long sample, PrintStream out) {
		mmm.calcForLogging(sample);
		
		
		
        double maxHeight=-1.0;
        for (int i = 0; i < mmm.treeIntervals.treeInput.get().size();i++) {
        	maxHeight = Math.max(mmm.treeIntervals.treeInput.get().get(i).getRoot().getHeight()+
        			mmm.treeIntervals.offset[i] + mmm.treeIntervals.rootLengthInput.get().get(i).getValue(), 
        			maxHeight);
        }
        String tree_string = "rem";

        for (int i = 0; i < mmm.treeIntervals.treeInput.get().size();i++) {
        	if (mmm.treeIntervals.treeInput.get().get(i).getExternalNodes().size()>=minClusterSizeInput.get()) {
        		Node root = getActualRoot(mmm.mappedTrees.get(i).getRoot());
        		String subtree_str = toNewick(root, root.getHeight());
        		subtree_str = subtree_str+";";
        		// remove the last bit of the tree string that indicates it
        		subtree_str = subtree_str.replace("]:0.0;", "");
        		subtree_str = subtree_str+ ",originHeight="+ (mmm.treeIntervals.rootLengthInput.get().get(i).getValue()+mmm.mappedTrees.get(i).getRoot().getHeight());

				subtree_str = subtree_str+ ",origin" + mmm.dynamics.typeTraitInput.getName() + "=" + mmm.dynamics.getStringStateValue((int) mmm.mappedTrees.get(i).getRoot().getLeft().getMetaData("location"));
        		subtree_str = subtree_str + "]:" + (maxHeight-root.getHeight()-mmm.treeIntervals.offset[i]);
        		tree_string = tree_string+"," +subtree_str;
        	}
        }
        tree_string = tree_string.replace("rem,", "(");
        tree_string = tree_string + "):0.0";
        tree_string = tree_string.replace("[&]", "");
        out.print("tree STATE_" + sample + " = ");
        out.print(tree_string);
        out.print(";");	        		
		
	}


	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
	
	private Node getActualRoot(Node node) {
		if (node.isLeaf() || node.getChildCount()==2) {
			return node;
		}else {
			return getActualRoot(node.getLeft()==null ? node.getRight() : node.getLeft());
		}
	}
	
	public String toNewick(Node node, double lastHeight) {
		StringBuffer buf = new StringBuffer();
		if (node.getLeft() != null) {
			if(node.getChildCount()==2) {
				buf.append("(");
				buf.append(toNewick(node.getLeft(), node.getHeight()));
				if (node.getRight() != null) {
					buf.append(',');
					buf.append(toNewick(node.getRight(), node.getHeight()));
				}
				buf.append(")");
			}else {
				buf.append(toNewick(node.getLeft(), lastHeight));
			}
		} else {
			if (!node.isLeaf())
				buf.append(node.getNr() + 1);
		}
		
		if(node.getChildCount()!=1) {

			if (!node.isLeaf()) {
				buf.append("[&obs=1,");
				buf.append(mmm.dynamics.typeTraitInput.getName() + "=" + mmm.dynamics.getStringStateValue((int) node.getMetaData("location")));			
			
				buf.append("]");
	
	
			} else {
	            buf.append(node.getID() + "[&obs=1,");
				
				if ( node.getMetaData("location")!=null) {
					buf.append(mmm.dynamics.typeTraitInput.getName() + "="
							+ mmm.dynamics.getStringStateValue((int) node.getMetaData("location")));
					buf.append("]");
				}
	
			}
	
			buf.append(":");
			mmm.appendDouble(buf, lastHeight-node.getHeight());
		}
		return buf.toString();
	}

}
