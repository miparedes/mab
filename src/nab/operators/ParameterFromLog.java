package nab.operators;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;


@Description("Scales a parameter or a complete beast.tree (depending on which of the two is specified.")
public class ParameterFromLog extends Operator {
    final public Input<List<StateNode>> upInput = new Input<>("up",
            "zero or more items to scale upwards", new ArrayList<>());

    public final Input<List<RealParameter>> parameterInput = new Input<>("parameter", "List of Parameters", new ArrayList<>());
    public final Input<String> logNameInput = new Input<>("logNames", "name of parameter in LogFile", Input.Validate.REQUIRED);
    public final Input<String> filenameInput = new Input<>("filename", "name of the log file", Input.Validate.REQUIRED);
    public final Input<Integer> burninInput = new Input<>("burnin", "burnin for the log file in percent", 10);


    
    Double[][] values; 

    @Override
    public void initAndValidate() {
    	
    	String[] useLabels = logNameInput.get().split("\\s+");
    	try {
			readLogFile(filenameInput.get(), (int) burninInput.get(), useLabels);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }



    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
    	int i=0;

        int counter = Randomizer.nextInt(values[0].length);

    	for (RealParameter param : parameterInput.get()) {
    		param.setValue(0, values[i][counter]);
    		i++;
    	}
    	counter++;
    	
    	return Double.POSITIVE_INFINITY;
    }
    
    @SuppressWarnings("unchecked")
	protected void readLogFile(String fileName, int burnInPercentage, String[] useLabels) throws IOException {
        System.out.println("\nLoading " + fileName);
        BufferedReader fin = new BufferedReader(new FileReader(fileName));
        String str;
        String m_sPreAmble = "";
        String[] m_sLabels = null;
        int data = 0;
        // first, sweep through the log file to determine size of the log
        while (fin.ready()) {
            str = fin.readLine();
            if (str.indexOf('#') < 0 && str.matches(".*[0-9a-zA-Z].*")) {
                if (m_sLabels == null)
                    m_sLabels = str.split("\\t");
                else
                    data++;
            } else {
                m_sPreAmble += str + "\n";
            }
        }
        int lines = Math.max(1, data / 80);
        int[] useItems = new int[useLabels.length];
        for (int i = 0; i < m_sLabels.length; i++) {
        	for (int j = 0; j < useLabels.length; j++) {
        		if (m_sLabels[i].contentEquals(useLabels[j])) {
        			useItems[j] = i;
        		}
        	}
        }
        
        int burnIn = data * burnInPercentage / 100;
        int total = data - burnIn;
        
        values = new Double[useItems.length][total];
        

        
        
        
        fin.close();
        fin = new BufferedReader(new FileReader(fileName));
        data = -burnIn - 1;
        // grab data from the log, ignoring burn in samples
        int reported = 0; 
        while (fin.ready()) {
            str = fin.readLine();
            if (str.indexOf('#') < 0 && str.matches("[-0-9].*")) {
                //data++;
                if (++data >= 0) {
                	String[] splitDat = str.split("\\s");       
                	for (int i = 0; i < values.length; i++)
                		values[i][data] = Double.parseDouble(splitDat[useItems[i]]);
                }
            }
        }
        fin.close();
    } // readLogFile



} // class ScaleOperator
