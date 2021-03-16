package nab.reassortment;



import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
public class AncestralFilteredSequenceExpandedLogger extends BEASTObject implements Function, Loggable {
	
    public Input<String> tagInput = new Input<String>("tag","label used to report trait", Validate.REQUIRED);

	final public Input<Boolean> logMutationsOnlyInput = new Input<>("logMutationsOnly", "if true, logs only mutations on edges instead of sequences on nodes", false);
	
	final public Input<List<AncestralStateExpandedTreeLikelihood>> ancestralTreeLikelihoodInput = new Input<>("ancestralTreeLikelihood", "ancestrsal tree likelihoods that," +
	"if more than one is inputted is assumed to be a filtered alignment", new ArrayList<>());
	
	final public Input<Boolean> translateInput = new Input<>("translate", "if true, dna is translated to amino acids", false);

	final public Input<Boolean> relativeToRootInput = new Input<>("relativeToRoot", "if true, difference is computed relative to root", false);

	
//	int [] siteStates;
	
    protected DataType dataType;
    
    int totalLength;
    List<Integer[]> mapping;
    
    Map<String, String> translationsTable;

	
	@Override
	public void init(PrintStream out) {
		ancestralTreeLikelihoodInput.get().get(0).updateTree();
		ancestralTreeLikelihoodInput.get().get(0).initialize();
		((Tree) ancestralTreeLikelihoodInput.get().get(0).tree).init(out);
		
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
				Integer[] newMap = new Integer[totalLength];
				for (int j = 0; j < newMap.length; j++)
					newMap[j] = j;
				mapping.add(newMap);
			}
		}
		if (translateInput.get())
			initTable();		
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


	boolean stop = false;
	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
			if (i==0) {
				ancestralTreeLikelihoodInput.get().get(i).updateTree();
			}else {
				ancestralTreeLikelihoodInput.get().get(i).setTree(ancestralTreeLikelihoodInput.get().get(0).tree);
			}
			ancestralTreeLikelihoodInput.get().get(i).calculateLogP();
		}
        out.print("tree STATE_" + sample + " = ");
        TreeInterface tree = ancestralTreeLikelihoodInput.get().get(0).tree;
        tree.getRoot().sort();

        out.print(toNewick(tree.getRoot(), null, null));
        out.print(";");
        if (stop) {
        	System.out.println(ancestralTreeLikelihoodInput.get().get(0).networkInput.get());
        	System.out.println(toNewick(tree.getRoot(), null, null));
        	System.exit(0);
        }
	}

	
    String toNewick(Node node, int[] parentSiteStates, List<String> translatedParent) {
        StringBuffer buf = new StringBuffer();
        
        int[] currSiteStates = new int[totalLength];
        List<String> translated = new ArrayList<>();
        
        for (int i = 0; i < currSiteStates.length; i++)
        	currSiteStates[i] =-1;
        
		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
		    int [] patternstates = ancestralTreeLikelihoodInput.get().get(i).getStatesForNode(
		    		ancestralTreeLikelihoodInput.get().get(0).tree, node);
		    
		    
		    for (int j = 0; j < mapping.get(i).length; j++) {
		    	currSiteStates[mapping.get(i)[j]] = patternstates[j];
		    }
		    
		    
		}
	    if (translateInput.get()) {
	    	translate(currSiteStates, translated);
	    	if (translatedParent==null)
	    		translatedParent = new ArrayList<>(translated);
	    }
	    if (relativeToRootInput.get()) {
	    	if (parentSiteStates==null) {
	    		parentSiteStates = new int[currSiteStates.length];
	    		System.arraycopy(currSiteStates, 0, parentSiteStates, 0, parentSiteStates.length);
	    	}
	    }

        if (node.getLeft() != null) {
            buf.append("(");
            if (relativeToRootInput.get())
            	buf.append(toNewick(node.getLeft(), parentSiteStates, translatedParent));
            else 
            	buf.append(toNewick(node.getLeft(), currSiteStates, translated));
            if (node.getRight() != null) {
                buf.append(',');
                if (relativeToRootInput.get())
                	buf.append(toNewick(node.getRight(), parentSiteStates, translatedParent));
                else
                	buf.append(toNewick(node.getRight(), currSiteStates, translated));
            }
            
            buf.append(")");
        } else {
            buf.append(node.getID());
        }
        
//		
//		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
//		    int [] patternstates = getStatesForNode(treeInput.get(), node);
//		    for (int j = 0; j < siteStates.length; j++) {
//		    	siteStates[j] = patternstates[dataInput.get().getPatternIndex(i)];
//		    }
		

	    if (logMutationsOnlyInput.get()) {
	    	if (parentSiteStates!=null) {
			    // look for differences
	    		if (!relativeToRootInput.get()) {
	    			if (!translateInput.get()) {
					    List<Integer> diffs = new ArrayList<>();
					    for (int i = 0; i < parentSiteStates.length; i++) {
					    	if (parentSiteStates[i]!=currSiteStates[i]) {
					    		diffs.add(i);
					    	}
					    	
					    }
					    if (diffs.size()>0) {
					    	buf.append("[&" + tagInput.get() + "=\"");
						    for (int i = 0; i < diffs.size(); i++) {
						    	buf.append(dataType.getCharacter(parentSiteStates[diffs.get(i)]) + "" + (diffs.get(i)+1) + 
						    			"" + dataType.getCharacter(currSiteStates[diffs.get(i)]));
						    	if (i<(diffs.size()-1)) {
						    		buf.append(",");
						    	}
						    }  
						    buf.append("\"");
					    	buf.append(",mutations="+diffs.size());
						    buf.append("]");
					    }
	    			}else {
					    List<Integer> diffs = new ArrayList<>();
					    for (int i = 0; i < translatedParent.size(); i++) {
					    	if (translatedParent.get(i)!=translated.get(i)) {
					    		diffs.add(i);
					    	}
					    	
					    }
					    if (diffs.size()>0) {
					    	buf.append("[&" + tagInput.get() + "=\"");
						    for (int i = 0; i < diffs.size(); i++) {
						    	buf.append(translatedParent.get(diffs.get(i)) + "" + (diffs.get(i)+1) + 
						    			"" + translated.get(diffs.get(i)));
						    	if (i<(diffs.size()-1)) {
						    		buf.append(",");
						    	}
						    }  
						    buf.append("\"");
					    	buf.append(",mutations="+diffs.size());
						    buf.append("]");
					    }

	    			}
	    		}else {
				    List<Integer> diffs = new ArrayList<>();
				    for (int i = 0; i < parentSiteStates.length; i++) {
				    	if (parentSiteStates[i]!=currSiteStates[i]) {
				    		diffs.add(i);
				    	}
				    	
				    }
				    if (diffs.size()>220) {
				    	System.out.println(Arrays.toString(parentSiteStates));				    	
				    	System.out.println(Arrays.toString(currSiteStates));
				    	System.out.println(node.getID());
				    	System.out.println(node.getNr());
				    	stop = true;
//				    	System.exit(0);
				    }
				    	
			    	buf.append("[&NT_muts=" + diffs.size());
	    			if (translateInput.get()) {
					    diffs = new ArrayList<>();
					    for (int i = 0; i < translatedParent.size(); i++) {
					    	if (translatedParent.get(i)!=translated.get(i)) {
					    		diffs.add(i);
					    	}
					    }
				    	buf.append(",AA_muts=" + diffs.size());
	    			}
				    buf.append("]");
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
        buf.append(Math.max(node.getLength(),0.000000000001));
        return buf.toString();
    }
    
	private void translate(int[] currSiteStates, List<String> translated) {
		for (int i = 0; i < currSiteStates.length-2; i=i+3) {
			String codon = dataType.getCharacter(currSiteStates[i]) +""+dataType.getCharacter(currSiteStates[i+1])+""+dataType.getCharacter(currSiteStates[i+2]); 
			translated.add(translationsTable.get(codon));
		}
//		System.out.println(translated);
//		System.exit(0);
		
	}


	@Override
	public void close(PrintStream out) {
		((Tree) ancestralTreeLikelihoodInput.get().get(0).tree).close(out);
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
	
	
	private void initTable() {
		String dna = "AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT";
		String aa = "K N K N T T T T R S R S I I M I Q H Q H P P P P R R R R L L L L E D E D A A A A G G G G V V V V O Y O Y S S S S O C W C L F L F";		
		
		translationsTable = new HashMap<>();
		String[] dnaarray = dna.split(" ");
		String[] aaarray = aa.split(" ");
		for (int i =0; i < dnaarray.length;i++)
			translationsTable.put(dnaarray[i], aaarray[i]);		
		
		
	}
}
