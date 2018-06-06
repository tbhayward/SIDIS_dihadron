/*
 * author Timothy B. Hayward
 * April 2018
 * SIDIS dihadron (dipion) event counter wtih various fitters
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

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*; // 

public class event_counter {
	// groovy program to calculate the rate of reconstructed dipion events following CLAS12 
	// detector simulation and reconstruction (gemc / myClara)
	// c.f. /volatile/clas12/thayward/SIDIS/dihadron/README

	private static double round (double value, int precision) {
		// just round a number to a certain precision
    int scale = (int) Math.pow(10, precision);
    return (double) Math.round(value * scale) / scale;
	}

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

		double beam_mass = 0.00051099894; // electron mass (GeV)
		double target_mass = 0.93827208; // proton mass (GeV)
		LorentzVector beam = new LorentzVector(0, 0, 
			Math.pow(Math.pow(10.6,2)-Math.pow(beam_mass,2),0.5), 10.6); // 10.6 GeV beam
		LorentzVector target = new LorentzVector(0, 0, 0, target_mass);

		GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6);
		GenericKinematicFitter rec_fitter = new event_builder_fitter(10.6);
		GenericKinematicFitter research_fitter = new enhanced_pid_fitter(10.6);
		EventFilter filter = new EventFilter("11:211:-211:X+:X-:Xn"); 
		// all events with electron, + and - pion and any number of other particles

		int event_counter = 0; // number of events analyzed
		int monte_carlo_counter = 0; // number of events with reconstructed MC particles
		int rec_counter = 0; // number of events reconstructed from generic CLAS12 fitter
		int research_counter = 0; // number of events reconstructed with research project fitter
		int channel_counter = 0; // number of events reconstructed with research project fitter and
			// meeting channel selection criteria

		for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files);
			// limit to a certain number of files defined by n_files
			HipoDataSource reader = new HipoDataSource();
			reader.open(hipo_list[current_file]); // open next hipo file

			while(reader.hasEvent()==true){ // cycle through events
				event_counter++; // new event analyzed
				HipoDataEvent event = reader.getNextEvent(); 
    			PhysicsEvent  mc_Event  = mc_fitter.getPhysicsEvent(event);
    			PhysicsEvent  rec_Event  = rec_fitter.getPhysicsEvent(event);
    			PhysicsEvent  research_Event  = research_fitter.getPhysicsEvent(event);

    			if(filter.isValid(mc_Event)==true){
    				monte_carlo_counter++; // monte carlo fitter returned recon dipion event
    			}
    			if(filter.isValid(rec_Event)==true){
    				rec_counter++; // generic CLAS12 fitter returned recon dipion event
    			}
    			if(filter.isValid(research_Event)==true){
    				research_counter++; // research fitter returned recon dipion event

    				Particle pipi = mc_Event.getParticle("[211]+[-211]");
    				Particle electron = mc_Event.getParticle("[11]");
    				if ((dis_test(beam, target, pipi, electron))&&
    					(sidis_test(beam, target, pipi, electron))) {
    					channel_counter++;
    				}
    			}
    		}
		}
		println(); println();
		println("Analyzed "+event_counter+" events.");
		println("monte_carlo_fitter reconstructed "+monte_carlo_counter+" events ("+
			round(100*monte_carlo_counter/event_counter,3)+"%)");
		println("CLAS12 Event Builder reconstructed "+rec_counter+" events ("+
			round(100*rec_counter/event_counter,3)+"%)");
		println("Enhanced PID fitter reconstructed "+research_counter+" events ("+
			round(100*research_counter/event_counter,3)+"%)");
		println("Enhanced PID fitter and channel filter reconstructed "+
			channel_counter+" events ("+round(100*channel_counter/event_counter,3)+"%)");
		println();
	}
}