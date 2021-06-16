/*
 * author Timothy B. Hayward
 * 
 * SIDIS dihadron 
 */

import java.io.File;

import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.*;
import org.jlab.clas12.physics.*;

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*; 
import analyzers.*;

// dilks CLAS QA analysis
import clasqa.QADB

// filetype for gathering files in directory
import groovy.io.FileType;


public class processing_single_hadrons {

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

    	def list = []

		def dir = new File(args[0])
		dir.eachFileRecurse (FileType.FILES) { file ->
		  list << file
		  // println(file.toString()); println(); 
		}

		println(); println(); println();

		String p1_Str;
		if (args.length < 2) {
			// assigns pi+ to p1
			println("WARNING: Specify a PDG PID for p1! Set to pi+. \n");
			p1_Str = "211";
		} else {
			p1_Str = args[1];
			println("Set p1 PID = "+p1_Str+"\n");
		}

		String output_file;
		if (args.length < 3) {
			// uses dummy name for output file if not specified
			println("WARNING: Specify an output file name. Set to \"dihadron_dummy_out.txt\".\n");
			output_file = "dihadron_dummy_out.txt"
		} else {
			output_file = args[2];
		}

		int n_files;
		if ((args.length < 4)||(Integer.parseInt(args[3])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			// n_files = hipo_list.size();
			n_files = list.size();
			println("There are "+hipo_list.size()+" or maybe "+list.size()+" number of files.")
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[2]);
		}

		File file = new File(output_file);
		file.bytes = new byte[0]

		int hadron_pair_counts = 0;
		GenericKinematicFitter research_fitter = new analysis_fitter(10.6041); // load my kinematic fitter/PID
		EventFilter filter = new EventFilter("11:"+p1_Str+":X+:X-:Xn"); // set filter for final states
		// setup QA database
		QADB qa = new QADB();

		int num_events = 0;
		int current_file = 0;
		// for (int current_file; current_file<n_files; current_file++) {
		while (current_file < n_files) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files); println(); println();
			// limit to a certain number of files defined by n_files

			HipoDataSource reader = new HipoDataSource();

			reader.open(list[current_file]); // open next hipo file
			current_file++;
			HipoDataEvent event = reader.getNextEvent(); 

			while(reader.hasEvent()==true){
				num_events++; 
				if (num_events%100000 == 0) { // not necessary
					print("processed: "+num_events+" events. ");
				}

				// get run and event numbers
				event = reader.getNextEvent();
			    int runnum = event.getBank("RUN::config").getInt('run',0); // collect info for QA
			    int evnum = event.getBank("RUN::config").getInt('event',0);

			    PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);

			    boolean process_event = false;
			    if (runnum == 11) { // if run number = 11 then it is MC and we don't use QA
			    	process_event = filter.isValid(research_Event);
			    } else {
			    	process_event = (filter.isValid(research_Event) && qa.OkForAsymmetry(runnum,evnum));
			    }
				if (process_event) {
					int num_p1 = research_Event.countByPid(p1_Str.toInteger());  // get # of particles w/ pid1

					for (int current_p1 = 0; current_p1 < num_p1; current_p1++) { // cycle over all combinations

						Hadron variables = new Hadron(event, research_Event, 
							p1_Str.toInteger(), current_p1);
						// this is my class for defining all relevant kinematic variables

						if (variables.channel_test(variables)) {
							int helicity = variables.get_helicity(); // helicity of event, might be 0

							// lab kinematics
							double e_p = variables.e_p();
							double e_theta = variables.e_theta();
							double e_phi = variables.e_phi();
							double p_p = variables.p_p();
							double p_theta = variables.p_theta();

							// DIS variables
							double Q2 = variables.Q2();
							double W = variables.W();
							double y = variables.y();
							double Mx = variables.Mx();

							// SIDIS variables
							double x = variables.x();
							double z = variables.z();
							double xF = variables.xF();
							double pT = variables.pT();
							double eta = variables.eta();
							double zeta = variables.zeta();

							// angles 
							double phi = variables.phi();

							// vertices 
							double vz_e = variables.vz_e();
							double vz_p = variables.vz_p();

							// append event to next line of the text file
							file.append(runnum+" "+evnum+" "+helicity+" ");
							file.append(e_p+" "+e_theta+" "+p_p+" "+p_theta+" ");
							file.append(Q2+" "+W+" "+Mx+" "+x+" "+y+" "+z+" "+xF+" ");
							file.append(pT+" "+zeta+" ");
							file.append(eta+" ");
							file.append(phi+" ");
							file.append(vz_e+" "+vz_p+" "+"\n");
						}
					}
				}
			}
			println(); println();
			print("1:runnum, 2:evnum, 3:helicity, 4:e_p, 5:e_theta, 6:p_p, 7:p_theta, ");
			print("8:Q2, 9:W, 10:Mx, 11:x, 12:y, 13:z, ");
			print("14:xF, ");
			print("15:pT, 16:zeta ");
			print("17:eta, ");
			print("18:phi, ");
			print("19: vz_e, 20: vz_p. \n");

			println(); println();
			println("Set p1 PID = "+p1_Str+"\n");
			println("output file is: "+file);
		}

	}
}