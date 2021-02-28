package nab.evolution.likelihood;



import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.evolution.datatype.DataType;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;

@Description("Reconstructs sequences at internal nodes and logs them in NEXUS format")
public class AncestralFilteredSequenceLogger extends BEASTObject implements Function, Loggable {
	
    public Input<String> tagInput = new Input<String>("tag","label used to report trait", Validate.REQUIRED);

	final public Input<Boolean> logMutationsOnlyInput = new Input<>("logMutationsOnly", "if true, logs only mutations on edges instead of sequences on nodes", false);
	
	final public Input<List<AncestralStateTreeLikelihood>> ancestralTreeLikelihoodInput = new Input<>("ancestralTreeLikelihood", "ancestrsal tree likelihoods that," +
	"if more than one is inputted is assumed to be a filtered alignment", new ArrayList<>());
	
//	int [] siteStates;
	
    protected DataType dataType;
    
    int totalLength;
    List<Integer[]> mapping;

	
	@Override
	public void init(PrintStream out) {
		
		
		((Tree) ancestralTreeLikelihoodInput.get().get(0).treeInput.get()).init(out);
		
		dataType =  ancestralTreeLikelihoodInput.get().get(0).dataInput.get().getDataType();
		totalLength = ancestralTreeLikelihoodInput.get().get(0).dataInput.get().getSiteCount();;

		mapping = new ArrayList<>();
		
		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
			if (ancestralTreeLikelihoodInput.get().get(i).dataInput.get() instanceof FilteredAlignment) {
				mapping.add(parseFilter(((FilteredAlignment) ancestralTreeLikelihoodInput.get().get(i).dataInput.get()).filterInput.get(),
						((FilteredAlignment) ancestralTreeLikelihoodInput.get().get(i).dataInput.get()).alignmentInput.get(),
						ancestralTreeLikelihoodInput.get().get(i).dataInput.get().getSiteCount()));
				
				totalLength = ((FilteredAlignment) ancestralTreeLikelihoodInput.get().get(i).dataInput.get()).alignmentInput.get().getSiteCount();
			}else {
				throw new IllegalArgumentException("only allows FilteredALignment based likelihood as input");
			}
		}
		
//		System.out.println(Arrays.toString(mapping.get(0)));
//		siteStates = new int[totalLength];	
	}
	
	
    private Integer[] parseFilter(String filterString, Alignment alignment, int redCount) {
        // parse filter specification
//        String filterString = filterInput.get();
        String[] filters = filterString.split(",");
        int[] from = new int[filters.length];
        int[] to = new int[filters.length];
        int[] step = new int[filters.length];
        for (int i = 0; i < filters.length; i++) {
            filterString = " " + filters[i] + " ";
            if (filterString.matches(".*-.*")) {
                // range, e.g. 1-100/3
                if (filterString.indexOf('\\') >= 0) {
                	String str2 = filterString.substring(filterString.indexOf('\\') + 1); 
                	step[i] = parseInt(str2, 1);
                	filterString = filterString.substring(0, filterString.indexOf('\\'));
                } else {
                	step[i] = 1;
                }
                String[] strs = filterString.split("-");
                from[i] = parseInt(strs[0], 1) - 1;
                to[i] = parseInt(strs[1], alignment.getSiteCount()) - 1;
            } else if (filterString.matches(".*:.*:.+")) {
                // iterator, e.g. 1:100:3
                String[] strs = filterString.split(":");
                from[i] = parseInt(strs[0], 1) - 1;
                to[i] = parseInt(strs[1], alignment.getSiteCount()) - 1;
                step[i] = parseInt(strs[2], 1);
            } else if (filterString.trim().matches("[0-9]*")) {
                from[i] = parseInt(filterString.trim(), 1) - 1;
                to[i] = from[i];
            	step[i] = 1;
            } else {
                throw new IllegalArgumentException("Don't know how to parse filter " + filterString);
            }
        }
        
//        System.out.println(Arrays.toString(from));
//        System.out.println(Arrays.toString(to));
//        System.out.println(Arrays.toString(step));
        
        
        boolean[] used = new boolean[alignment.getSiteCount()];
        for (int i = 0; i < to.length; i++) {
            for (int k = from[i]; k <= to[i]; k += step[i]) {
            	used[k] = true;
            }
        }
//        System.out.println(Arrays.toString(map));
//        System.exit(0);
        
        Integer[] map = new Integer[redCount];
        int c=0;
        for (int i = 0; i < used.length; i++) {
        	if (used[i]) {
        		map[c] = i;
        		c++;
        	}
        }

        return map;
    }
    
    int parseInt(String str, int defaultValue) {
        str = str.replaceAll("\\s+", "");
        try {
            return Integer.parseInt(str);
        } catch (Exception e) {
            return defaultValue;
        }
    }


	
	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
			ancestralTreeLikelihoodInput.get().get(i).calculateLogP();
			ancestralTreeLikelihoodInput.get().get(i).redrawAncestralStates();

		}
        out.print("tree STATE_" + sample + " = ");
        TreeInterface tree = ancestralTreeLikelihoodInput.get().get(0).treeInput.get();
        tree.getRoot().sort();
        out.print(toNewick(tree.getRoot(), null));
        out.print(";");
	}

	
    String toNewick(Node node, int[] parentSiteStates) {
        StringBuffer buf = new StringBuffer();
        
        int[] currSiteStates = new int[totalLength];
        
        for (int i = 0; i < currSiteStates.length; i++)
        	currSiteStates[i] =-1;
        
		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
		    int [] patternstates = ancestralTreeLikelihoodInput.get().get(i).getStatesForNode(
		    		ancestralTreeLikelihoodInput.get().get(0).treeInput.get(), node);
		    
		    
		    for (int j = 0; j < mapping.get(i).length; j++) {
		    	currSiteStates[mapping.get(i)[j]] = patternstates[ancestralTreeLikelihoodInput.get().get(i).dataInput.get().getPatternIndex(j)];
		    }
		}

        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), currSiteStates));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), currSiteStates));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }
        
		
//		System.out.println(ancestralTreeLikelihoodInput.get().get(0).dataInput.get().getPatternIndex(578));
//		System.out.println(Arrays.toString(ancestralTreeLikelihoodInput.get().get(0).dataInput.get().getPattern(0)));
//		System.out.println(Arrays.toString(ancestralTreeLikelihoodInput.get().get(0).dataInput.get().getPattern(1)));
//		System.out.println(Arrays.toString(ancestralTreeLikelihoodInput.get().get(0).dataInput.get().getPattern(578)));
//		
//		System.out.println(Arrays.toString(mapping.get(0)));
//		System.out.println(Arrays.toString(siteStates));
//
//		System.out.println(siteStates);
//		
//		System.exit(0);
//		
//		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
//		    int [] patternstates = getStatesForNode(treeInput.get(), node);
//		    for (int j = 0; j < siteStates.length; j++) {
//		    	siteStates[j] = patternstates[dataInput.get().getPatternIndex(i)];
//		    }
		

	    if (logMutationsOnlyInput.get()) {
	    	if (parentSiteStates!=null) {
			    // look for differences
			    List<Integer> diffs = new ArrayList<>();
			    for (int i = 0; i < parentSiteStates.length; i++) {
			    	if (parentSiteStates[i]!=currSiteStates[i]) {
			    		diffs.add(i);
			    	}
			    	
			    }
			    if (diffs.size()>0) {
			    	buf.append("[&" + tagInput.get() + "=\"");
				    for (int i = 0; i < diffs.size(); i++) {
				    	buf.append(dataType.getCharacter(parentSiteStates[diffs.get(i)]) + ""+diffs.get(i)+ 
				    			"" + dataType.getCharacter(currSiteStates[diffs.get(i)]));
				    	if (i<(diffs.size()-1)) {
				    		buf.append(",");
				    	}
				    }
				    buf.append("\"]");
			    }
	
	//	    	buf.append("[&");
	//	    	for (int k = 0; k < seq.length(); k++) {
	//	    		buf.append((k > 0 ? "," : "") + tagInput.get()
	//	    		+ (k < 10 ? "0":"")
	//	    		+ (k < 100 ? "0":"")
	//	    		+ k + "=\"" + seq.charAt(k) + "\"");
	//	    	}
	    	}
	    	
	    } else {
		    String seq_current = dataType.encodingToString(currSiteStates);
	    	buf.append("[&" + tagInput.get() + "=\"" + seq_current + "\"]");
	    }

	    buf.append(':');
        buf.append(node.getLength());
        return buf.toString();
    }
    
	@Override
	public void close(PrintStream out) {
		((Tree)ancestralTreeLikelihoodInput.get().get(0).treeInput.get()).close(out);
	}

	@Override
	public int getDimension() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getArrayValue(int dim) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}	

}
