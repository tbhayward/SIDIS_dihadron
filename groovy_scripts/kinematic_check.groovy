/*
 * author Timothy B. Hayward
 * May 2018
 * SIDIS dihadron (dipion) kinematic check
 * just check basic kinematic distributions
 */

import java.io.File;

import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.*;
import org.jlab.clas12.physics.*;

import javax.swing.JFrame;
import org.jlab.groot.data.*;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.io.base.DataEvent;

import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.fitter.DataFitter;

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*; 
import histogram_calculations.*;

public class resolution {

	public static double[] bins(int n_bins, double min_bin, double max_bin) {
		double[] bin_centers = new double[n_bins-1];
		double running_center = min_bin+(0.5)*(max_bin-min_bin)/n_bins;
		for (int i=0; i<n_bins-1; i++) {
			bin_centers[i] = running_center;
			running_center+= (max_bin-min_bin)/n_bins;
		}

		return bin_centers;
	}

	public static F1D fit_function(int n_bins, double min_bin, double max_bin, H1F histogram,
		double mu, double sigma) {

		// first fit to a Gaussian and check chi^2/dof 
		String func_string = "[amp]*gaus(x,[mean],[sigma])";
		F1D func = new F1D("gauss",func_string,min_bin,max_bin);
		func.setParameter(0, 1.0);
		func.setParameter(1, mu);
		func.setParameter(2, sigma);
		DataFitter.fit(func, histogram, "Q"); //No options uses error for sigma
		double chi2dof = func.getChiSquare()/(n_bins-1);

		// add polynomials until chi2dof stops improving (<1 or >chi2dof)
		boolean chi2dof_check = true; 
		F1D test_func = new F1D("gauss+?",func_string,min_bin,max_bin);
		test_func.setParameter(0, 1.0);
		test_func.setParameter(1, mu);
		test_func.setParameter(2, sigma);
		int poly_order = 0; // current order of polynomial added to func
		while(chi2dof_check) {
		// for (int current_order=0; current_order < 3; current_order++) {
			func_string+="+[p"+poly_order+"]"
			for (int i=0; i<poly_order; i++) {
				func_string+="*x";
				test_func = new F1D("gauss+?",func_string,min_bin,max_bin);
				test_func.setParameter(i+3, 1.0);
			}
			test_func.setParameter(0, 1.0);
			test_func.setParameter(1, mu);
			test_func.setParameter(2, sigma);
			DataFitter.fit(test_func, histogram, "Q"); //No options uses error for sigma

			if ((test_func.getChiSquare()/(n_bins-1))<1) {
				println("Chi2/dof dropped below 1.0. Ending fitting routine.");
				chi2dof_check = false;
			} else if ((test_func.getChiSquare()/(n_bins-1))>chi2dof) {
				println("Chi2/dof of gaus+p(i+1) > gaus+p(i). Ending fitting routine.");
				chi2dof_check = false;
			} else {
				chi2dof = test_func.getChiSquare()/(n_bins-1);
				func = test_func;
			}
			poly_order++;	
		}
		println(chi2dof);
		return func;
	}

	public static void main(String[] args) {
		File[] hipo_list;
		if (args.length == 0) {
			// exits program if input directory not specified 
        	println("ERROR: Please enter a hipo file directory as the first argument");
       		System.exit(0);
    	} else {
    		File directory = new File(args[0]);
    		hipo_list = directory.listFiles();
    	}

    	int n_files;
		if ((args.length < 2)||(Integer.parseInt(args[1])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			n_files = hipo_list.size();
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[1]);
		}


		int n_bins = 200;
		double min_bin = -1;
		double max_bin = 10;

		// JFrame frame = new JFrame("Timothy's u-u TMDGen Results");
		// JFrame frame = new JFrame("Harut's l-u Pythia Results");
		// JFrame frame = new JFrame("Timothy's MC u-u TMDGen Results");
		JFrame frame = new JFrame("");
		frame.setSize(900,470);

		double[] bin_centers = new double[n_bins-1];
		bin_centers = bins(n_bins, min_bin, max_bin);
		EmbeddedCanvas canvas = new EmbeddedCanvas();	
		H1F histogram = new H1F("",n_bins,min_bin,max_bin);
		histogram.setTitleX("pid #chi^2");
		histogram.setTitleY("Counts (Normalized)");

		GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6);
		GenericKinematicFitter rec_fitter = new event_builder_fitter(10.6);
		GenericKinematicFitter research_fitter = new enhanced_pid_fitter(10.6);
		EventFilter filter = new EventFilter("11:X+:X-:Xn"); 
		// all events with electron, + and - pion and any number of other particles

		for (int current_file; current_file<n_files; current_file++) {
		// for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files);
			// limit to a certain number of files defined by n_files
			HipoDataSource reader = new HipoDataSource();
			reader.open(hipo_list[current_file]); // open next hipo file

			while(reader.hasEvent()==true){ // cycle through events
				HipoDataEvent event = reader.getNextEvent(); 
    			PhysicsEvent  mc_Event  = mc_fitter.getPhysicsEvent(event);
    			PhysicsEvent  rec_Event  = rec_fitter.getPhysicsEvent(event);
    			PhysicsEvent  research_Event  = research_fitter.getPhysicsEvent(event);

    			double beam_mass = 0.00051099894; // electron mass (GeV)
				double target_mass = 0.93827208; // proton mass (GeV)
				LorentzVector beam = new LorentzVector(0, 0, 
					Math.pow(Math.pow(10.6,2)-Math.pow(beam_mass,2),0.5), 10.6); // 10.6 GeV beam
				LorentzVector target = new LorentzVector(0, 0, 0, target_mass)

    			if (filter.isValid(research_Event)) {
    	// 			Particle electron = research_Event.getParticle("[11]");
    	// 			Particle pipi = research_Event.getParticle("[211]+[-211]");

					// double Q2 = 4*beam.e()*electron.e()*Math.pow(electron.theta()/2,2);
					// double first_term = target.e()*target.e();
					// double second_term = 2*target.e()*(beam.e()-electron.e());
					// double third_term = - 4*beam.e()*electron.e()*Math.pow(Math.sin(electron.theta()/2),2);
					// double W = Math.pow(first_term + second_term + third_term, 0.5);

					// LorentzVector missing_part = new LorentzVector(pipi.px()+electron.px(), 
    	// 				pipi.py()+electron.py(), beam.pz()-pipi.pz()-electron.pz(), 
    	// 				beam.e()+target.e()-pipi.e()-electron.e());
    	// 			double missing_mass = Math.pow(Math.pow(missing_part.e(),2)-
    	// 				Math.pow(missing_part.p(),2),0.5);
    				HipoDataBank recBank= (HipoDataBank) event.getBank("REC::Particle");
    				for (int i=0; i<recBank.rows(); i++){
    					int pid = recBank.getInt("pid", i);
    					if (pid==11) {
    						float chi2pid = recBank.getFloat("chi2pid", i);
    						histogram.fill(chi2pid);
    					}
    				}
    			}
			}
		}

		histogram.unit();
		// draw the histogram
		canvas.setAxisLabelSize(14);
		canvas.setStatBoxFontSize(14);
		canvas.setAxisTitleSize(18);
		canvas.setTitleSize(18);
		canvas.draw(histogram);
		frame.add(canvas);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		canvas.setStatBoxFontSize(18);
		histogram.setOptStat(10000);
    }
}