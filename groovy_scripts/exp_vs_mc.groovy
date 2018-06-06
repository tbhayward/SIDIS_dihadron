/*
 * author Timothy B. Hayward
 * May 2018
 * SIDIS dihadron (dipion) REC::Particle vs MC::Particle
 * just check basic kinematic distributions and compare
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

	public static boolean dis_test(LorentzVector beam, LorentzVector target, 
			Particle pipi, Particle electron) {
    	double Q2 = 4*beam.e()*electron.e()*Math.pow(electron.theta()/2,2);

    	// target.e() is target mass in this case
    	double first_term = target.e()*target.e();
		double second_term = 2*target.e()*(beam.e()-electron.e());
		double third_term = - 4*beam.e()*electron.e()*Math.pow(Math.sin(electron.theta()/2),2);
		double W = Math.pow(first_term + second_term + third_term, 0.5);

		if ((Q2>1.0)&(W>2.0)) {
			return true;
		} else {
			return false;
		}
	}

	public static boolean sidis_test(LorentzVector beam, LorentzVector target, 
			Particle pipi, Particle electron) {
    	LorentzVector missing_part = new LorentzVector(pipi.px()+electron.px(), 
    		pipi.py()+electron.py(), beam.pz()-pipi.pz()-electron.pz(), 
    		beam.e()+target.e()-pipi.e()-electron.e());
    	double missing_mass = Math.pow(Math.pow(missing_part.e(),2)-
    		Math.pow(missing_part.p(),2),0.5);

		if (missing_mass>1.05) {
			return true;
		} else {
			return false;
		}
	}

	public static double[] bins(int n_bins, double min_bin, double max_bin) {
		double[] bin_centers = new double[n_bins-1];
		double running_center = min_bin+(0.5)*(max_bin-min_bin)/n_bins;
		for (int i=0; i<n_bins-1; i++) {
			bin_centers[i] = running_center;
			running_center+= (max_bin-min_bin)/n_bins;
		}

		return bin_centers;
	}

	public static double[] histogram_values(int n_bins, H1F histogram) {
		double[] histogram_values = new double[n_bins-1];
		for (int i=0; i<n_bins-1; i++) {
			histogram_values[i] = histogram.getBinContent(i);
		}

		return histogram_values;
	}

	public static int normalization_constant(int n_bins, H1F histogram) {
		int runSum=0;
		for (int i=0; i<n_bins-1; i++) {
			runSum+=histogram.getBinContent(i);
		}
		return runSum;
	}

	public static void main(String[] args) {
		File[] hipo_list_mc;
		if (args.length == 0) {
			// exits program if input directory not specified 
        	println("ERROR: Please enter a hipo file directory as the first argument");
       		System.exit(0);
    	} else {
    		File directory = new File(args[0]);
    		hipo_list_mc = directory.listFiles();
    	}

    	File[] hipo_list_exp;
		if (args.length < 2) {
			// exits program if input directory not specified 
        	println("ERROR: Please enter a hipo file directory as the second argument");
       		System.exit(1);
    	} else {
    		File directory = new File(args[1]);
    		hipo_list_exp = directory.listFiles();
    	}

    	int n_files_mc;
		if ((args.length < 3)||(Integer.parseInt(args[2])>hipo_list_mc.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			n_files_mc = hipo_list_mc.size();
		} else{
			// if specified, convert to int
			n_files_mc = Integer.parseInt(args[2]);
		}

		int n_files_exp;
		if ((args.length < 4)||(Integer.parseInt(args[3])>hipo_list_exp.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			n_files_exp = hipo_list_exp.size();
		} else{
			// if specified, convert to int
			n_files_exp = Integer.parseInt(args[3]);
		}


		int n_bins = 100;
		double min_bin = 0;
		double max_bin = 2;

		JFrame frame = new JFrame("");
		frame.setSize(900,470);

		double[] bin_centers = new double[n_bins-1];
		bin_centers = bins(n_bins, min_bin, max_bin);
		EmbeddedCanvas canvas = new EmbeddedCanvas();	
		H1F rec_histogram = new H1F("",n_bins,min_bin,max_bin);
		H1F mc_histogram = new H1F("",n_bins,min_bin,max_bin);
		// histogram.setTitleX("Electron Energy (GeV)");
		// histogram.setTitleY("Counts (Normalized)");

		GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6);
		GenericKinematicFitter rec_fitter = new event_builder_fitter(10.6);
		GenericKinematicFitter research_fitter = new enhanced_pid_fitter(10.6);
		EventFilter filter = new EventFilter("11:211:-211:X+:X-:Xn");
		// all events with electron, + and - pion and any number of other particles

		double beam_mass = 0.00051099894; // electron mass (GeV)
		double target_mass = 0.93827208; // proton mass (GeV)
		LorentzVector beam = new LorentzVector(0, 0, 
			Math.pow(Math.pow(10.6,2)-Math.pow(beam_mass,2),0.5), 10.6); // 10.6 GeV beam
		LorentzVector target = new LorentzVector(0, 0, 0, target_mass);

		for (int current_file; current_file<n_files_mc; current_file++) {
		// for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files_mc+" monte carlo files");
			// limit to a certain number of files defined by n_files
			HipoDataSource reader = new HipoDataSource();
			reader.open(hipo_list_mc[current_file]); // open next hipo file

			while(reader.hasEvent()==true){ // cycle through events
				HipoDataEvent event = reader.getNextEvent(); 
    			PhysicsEvent  mc_Event  = mc_fitter.getPhysicsEvent(event);

    			if(filter.isValid(mc_Event)==true){ 

    				Particle pipi = mc_Event.getParticle("[211]+[-211]");
    				Particle electron = mc_Event.getParticle("[11]");
    				if ((dis_test(beam, target, pipi, electron))&
    					(sidis_test(beam, target, pipi, electron)) ) {
    					double Q2 = 4*10.6*electron.e()*Math.pow(electron.theta()/2,2);
    					double xb = Q2 / (2*target.e()*(beam.e()-electron.e()));
    					double z = pipi.e()/(beam.e()-electron.e());
    					double mpipi = pipi.mass();
    					mc_histogram.fill(mpipi);
    				}
    			}
			}
		}

		for (int current_file; current_file<n_files_exp; current_file++) {
		// for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files_exp+" experimental files");
			// limit to a certain number of files defined by n_files
			HipoDataSource reader = new HipoDataSource();
			reader.open(hipo_list_exp[current_file]); // open next hipo file

			while(reader.hasEvent()==true){ // cycle through events
				HipoDataEvent event = reader.getNextEvent(); 
    			PhysicsEvent  rec_Event  = rec_fitter.getPhysicsEvent(event);
    			PhysicsEvent  research_Event  = research_fitter.getPhysicsEvent(event);

    			if(filter.isValid(research_Event)==true){ 
    				Particle pipi = research_Event.getParticle("[211]+[-211]");
    				Particle electron = research_Event.getParticle("[11]");
    				if ((dis_test(beam, target, pipi, electron))&
    					(sidis_test(beam, target, pipi, electron)) ) {
    					double Q2 = 4*10.6*electron.e()*Math.pow(electron.theta()/2,2);
    					double xb = Q2 / (2*target.e()*(beam.e()-electron.e()));
    					double z = pipi.e()/(beam.e()-electron.e());
    					double mpipi = pipi.mass();
    					rec_histogram.fill(mpipi);
    				}
    			}
			}
		}

		double[] rec_values = histogram_values(n_bins, rec_histogram);
		double[] mc_values = histogram_values(n_bins, mc_histogram);
		int rec_runSum = normalization_constant(n_bins, rec_histogram);
		int mc_runSum = normalization_constant(n_bins, mc_histogram);
		GraphErrors rec_points = new GraphErrors();
		GraphErrors mc_points = new GraphErrors();
		for (int i=0; i<n_bins-1; i++) {
			rec_points.addPoint(bin_centers[i], rec_values[i]/rec_runSum, 0, 
				Math.pow(rec_runSum,-1)*Math.pow(rec_values[i],0.5));
			mc_points.addPoint(bin_centers[i], mc_values[i]/mc_runSum, 0, Math.pow(mc_runSum,-1)*
				Math.pow(mc_values[i],0.5));
		}

		rec_points.setMarkerColor(2);
		rec_points.setLineColor(2);
		rec_points.setMarkerStyle(0);
		rec_points.setMarkerSize(4);
		rec_points.setLineThickness(2);

		mc_points.setMarkerColor(4);
		mc_points.setLineColor(4);
		mc_points.setMarkerStyle(0);
		mc_points.setMarkerSize(4);
		mc_points.setLineThickness(2);
		
		
		rec_points.setTitleX("xB");
		rec_points.setTitleY("Normalized Counts");
		canvas.getPad(0).setTitle("");
		canvas.setTitleSize(38);
		canvas.setAxisTitleSize(26);
		canvas.setAxisLabelSize(22);
		canvas.setStatBoxFontSize(18);
		canvas.draw(rec_points);
		canvas.draw(mc_points,"same");


		frame.add(canvas);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);



		// draw the histogram
		// canvas.draw(histogram);
		// frame.add(canvas);
		// frame.setLocationRelativeTo(null);
		// frame.setVisible(true);
		// canvas.setStatBoxFontSize(18);
		// histogram.setOptStat(1111);
    }
}