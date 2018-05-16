/*
 * author Timothy B. Hayward
 * May 2018
 * SIDIS dihadron (dipion) resolution tests
 * investigates difference between MC::Particle and REC::Particle banks
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
		double min_bin = -0.2;
		double max_bin = 0.2;

		JFrame frame = new JFrame("some frame here");
		frame.setSize(900,470);

		double[] bin_centers = new double[n_bins-1];
		bin_centers = bins(n_bins, min_bin, max_bin);
		EmbeddedCanvas canvas = new EmbeddedCanvas();	
		H1F histogram = new H1F("Resolution",n_bins,min_bin,max_bin);
		histogram.setTitleX("MC Electron Energy - REC Electron Energy (Gev)");
		histogram.setTitleY("Counts (Normalized)");

		GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6);
		GenericKinematicFitter generic_fitter = new generic_rec_fitter(10.6);
		// GenericKinematicFitter research_fitter = new SIDIS_fitter(10.6);
		EventFilter filter = new EventFilter("11:211:-211:X+:X-:Xn"); 
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
    			PhysicsEvent  generic_Event  = generic_fitter.getPhysicsEvent(event);
    			// PhysicsEvent  research_Event  = research_fitter.getPhysicsEvent(event);

    			if(filter.isValid(generic_Event)==true){ 
    			// only interested if REC::Particle has reconstructed event
    			Particle mc_e = mc_Event.getParticle("[11]");
    			Particle rec_e = generic_Event.getParticle("[11]");
    			histogram.fill(mc_e.e()-rec_e.e());

    			}
			}
		}

		histogram.unit();
		double mu = bin_centers[histogram.getMaximumBin()];
		int half_height_index = 0;
		double half_height_test = 0;
		while(half_height_test < 0.5) {
			half_height_test = histogram.getDataY(half_height_index);
			half_height_index++;
		}
		double sigma = mu-bin_centers[half_height_index];
		// println(mu+" +/- "+sigma);

		F1D func = fit_function(n_bins, min_bin, max_bin, histogram, mu, sigma);
		// draw the histogram
		canvas.draw(histogram);
		frame.add(canvas);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		canvas.setStatBoxFontSize(18);
		histogram.setOptStat(10);
		func.show();
		func.setLineColor(2); func.setLineWidth(5); func.setLineStyle(0);
		canvas.draw(func,"same");
		func.setOptStat(1110);
		println("Fit has chi2dof = "+(func.getChiSquare()/(n_bins-1)));
    }
}