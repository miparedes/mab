package nab.util;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Tree;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

@Description("Logs tree heights and origins using trees and offset values")

public class MultiTreeStatLogger extends CalculationNode implements Loggable {
    final public Input<Tree> treeInput = new Input<>("tree", "tree for which to calculate the intervals", Input.Validate.REQUIRED);
    final public Input<RealParameter> offsetInput = new Input<>("offset", "time offset ", Input.Validate.REQUIRED);
    final public Input<RealParameter> rootLengthInput = new Input<>("rootLength", "time offset ", Input.Validate.REQUIRED);

    final public Input<Boolean> heightOnlyInput = new Input<>("heightOnly", "time offset ", false);
    final public Input<Boolean> originOnlyInput = new Input<>("originOnly", "time offset ", false);

		
	int nrIntervals,dims;
	boolean NeTimeVariant, migTimeVariant;
	
	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void init(PrintStream out) {
		if (heightOnlyInput.get()) {
			out.print(String.format("%s.trueHeight\t", treeInput.get().getID()));
		}else if (originOnlyInput.get()) {
			out.print(String.format("%s.trueOrigin\t", treeInput.get().getID()));
		}else {
			out.print(String.format("%s.trueHeight\t", treeInput.get().getID()));
			out.print(String.format("%s.trueOrigin\t", treeInput.get().getID()));
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		if (heightOnlyInput.get()) {
			out.print((treeInput.get().getRoot().getHeight() + offsetInput.get().getValue())  + "\t");
		}else if (originOnlyInput.get()) {
			out.print((treeInput.get().getRoot().getHeight() + rootLengthInput.get().getValue() + offsetInput.get().getValue())  + "\t");
		}else {
			out.print((treeInput.get().getRoot().getHeight() + offsetInput.get().getValue())  + "\t");
			out.print((treeInput.get().getRoot().getHeight() + rootLengthInput.get().getValue() + offsetInput.get().getValue())  + "\t");
		}
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}



}
