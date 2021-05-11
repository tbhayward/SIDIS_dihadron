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

public class processing_mc {

	public static double phi_calculation (double x, double y) {
		// tracks are given with Cartesian values and so must be converted to cylindrical
		double phi = Math.toDegrees(Math.atan2(x,y));
		phi = phi - 90;
		if (phi < 0) {
			phi = 360 + phi;
		}
		phi = 360 - phi;
		return phi;	
	}

	public static double theta_calculation (double x, double y, double z) {
		// convert cartesian coordinates to polar angle
		double r = Math.pow(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2),0.5);
		return (double) (180/Math.PI)*Math.acos(z/r);
	}



	public static void main(String[] args) {

		double scale = 2.0;


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

		// list.each {
		// 	println it.path
		// }

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

		String p2_Str;
		if (args.length < 3) {
			// assigns pi- to p2
			println("WARNING: Specify a PDG PID for p2! Set to pi-. \n");
			p2_Str = "-211";
		} else {
			p2_Str = args[2];
			println("Set p2 PID = "+p2_Str+"\n");
		}

		String output_file;
		if (args.length < 4) {
			// uses dummy name for output file if not specified
			println("WARNING: Specify an output file name. Set to \"dihadron_dummy_out.txt\".\n");
			output_file = "dihadron_dummy_out.txt"
		} else {
			output_file = args[3];
		}

		int n_files;
		if ((args.length < 5)||(Integer.parseInt(args[4])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			// n_files = hipo_list.size();
			n_files = list.size();
			println("There are "+hipo_list.size()+" or maybe "+list.size()+" number of files.")
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[4]);
		}

		File file = new File(output_file);
		file.bytes = new byte[0]

		int hadron_pair_counts = 0;
		GenericKinematicFitter research_fitter = new analysis_fitter(10.6041);
		GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6041);
		EventFilter filter = new EventFilter("11:"+p1_Str+":"+p2_Str+":X+:X-:Xn");
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
				if (num_events%100000 == 0) { 
					print("processed: "+num_events+" events.     ");
				}

				// get run and event numbers
				event = reader.getNextEvent();
			    int runnum = event.getBank("RUN::config").getInt('run',0);
			    int evnum = event.getBank("RUN::config").getInt('event',0);

			    PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);
			    PhysicsEvent mc_Event = mc_fitter.getPhysicsEvent(event);

			    boolean process_event = false;
			    if (runnum == 11) {
			    	process_event = filter.isValid(research_Event);
			    } else {
			    	println("runnum neq 11, is this MC? Script is only designed for MC.")
			    	process_event = (filter.isValid(research_Event) && qa.OkForAsymmetry(runnum,evnum));
			    }
				if (process_event) {

					HipoDataBank recBank = (HipoDataBank) event.getBank("REC::Event");
					HipoDataBank lundBank = (HipoDataBank) event.getBank("MC::Lund");
					HipoDataBank mcBank = (HipoDataBank) event.getBank("MC::Particle");

					int num_p1 = research_Event.countByPid(p1_Str.toInteger()); 
					int num_p2 = research_Event.countByPid(p2_Str.toInteger());

					for (int current_p1 = 0; current_p1 < num_p1; current_p1++) {
						for (int current_p2 = 0; current_p2 < num_p2; current_p2++) {

							Particle exp_e = research_Event.getParticleByPid(11,0);
							Particle exp_p1 = research_Event.getParticleByPid(p1_Str.toInteger(),current_p1);
							Particle exp_p2 = research_Event.getParticleByPid(p2_Str.toInteger(),current_p2);

							Dihadrons variables = new Dihadrons(event, research_Event, 
								p1_Str.toInteger(), current_p1, p2_Str.toInteger(), current_p2);
							Dihadrons mc_variables = new Dihadrons(event, mc_Event, 
								p1_Str.toInteger(), current_p1, p2_Str.toInteger(), current_p2);

							if (variables.channel_test(variables)) {

								// lab kinematics data
								double e_p = variables.e_p();
								double e_theta = variables.e_theta();
								double p1_p = variables.p1_p();
								double p1_theta = variables.p1_theta();
								double p2_p = variables.p2_p();
								double p2_theta = variables.p2_theta();

								// lab kinematics MC
								double mc_e_p = mc_variables.e_p();
								double mc_e_theta = mc_variables.e_theta();
								double mc_p1_p = mc_variables.p1_p();
								double mc_p1_theta = mc_variables.p1_theta();
								double mc_p2_p = mc_variables.p2_p();
								double mc_p2_theta = mc_variables.p2_theta();

								// DIS variables data
								double Q2 = variables.Q2();
								double W = variables.W();
								double y = variables.y();
								double Mx = variables.Mx();
								double Mx1 = variables.Mx1();
								double Mx2 = variables.Mx2();

								// DIS variables MC
								double mc_Q2 = mc_variables.Q2();
								double mc_W = mc_variables.W();
								double mc_y = mc_variables.y();
								double mc_Mx = mc_variables.Mx();
								double mc_Mx1 = mc_variables.Mx1();
								double mc_Mx2 = mc_variables.Mx2();

								// SIDIS variables data
								double x = variables.x();
								double z = variables.z();
								double xF = variables.xF();
								double pT = variables.pT();
								double eta = variables.eta();
								double eta_gN = variables.eta_gN();

								// SIDIS variables MC
								double mc_x = mc_variables.x();
								double mc_z = mc_variables.z();
								double mc_xF = mc_variables.xF();
								double mc_pT = mc_variables.pT();
								double mc_eta = mc_variables.eta();
								double mc_eta_gN = mc_variables.eta_gN();

								// SIDIS dihadron variables data
								double z1 = variables.z1();
								double z2 = variables.z2();
								double xF1 = variables.xF1();
								double xF2 = variables.xF2();
								double zeta = variables.zeta(); // assumption here is p1 is a proton
								// otherwise zeta is meaningless
								double Mh = variables.Mh();
								double pT1 = variables.pT1();
								double pT2 = variables.pT2();
								double pTpT = variables.pTpT();
								double eta1 = variables.eta1();
								double eta2 = variables.eta2();
								double Delta_eta = variables.Delta_eta();
								double eta1_gN = variables.eta1_gN();
								double eta2_gN = variables.eta2_gN();

								// SIDIS dihadron variables MC
								double mc_z1 = mc_variables.z1();
								double mc_z2 = mc_variables.z2();
								double mc_xF1 = mc_variables.xF1();
								double mc_xF2 = mc_variables.xF2();
								double mc_zeta = mc_variables.zeta(); // assumption here is p1 is a proton
								// otherwise zeta is meaningless
								double mc_Mh = mc_variables.Mh();
								double mc_pT1 = mc_variables.pT1();
								double mc_pT2 = mc_variables.pT2();
								double mc_pTpT = mc_variables.pTpT();
								double mc_eta1 = mc_variables.eta1();
								double mc_eta2 = mc_variables.eta2();
								double mc_Delta_eta = mc_variables.Delta_eta();
								double mc_eta1_gN = mc_variables.eta1_gN();
								double mc_eta2_gN = mc_variables.eta2_gN();

								// angles data
								double phi1 = variables.phi1();
								double phi2 = variables.phi2();
								double Delta_phi = variables.Delta_phi();
								double phih = variables.phih();
								double phiR = variables.phiR();
								double theta = variables.theta();

								// angles MC
								double mc_phi1 = mc_variables.phi1();
								double mc_phi2 = mc_variables.phi2();
								double mc_Delta_phi = mc_variables.Delta_phi();
								double mc_phih = mc_variables.phih();
								double mc_phiR = mc_variables.phiR();
								double mc_theta = mc_variables.theta();

								// depolarization factors data
								double Depolarization_A = variables.Depolarization_A();
								double Depolarization_B = variables.Depolarization_B();
								double Depolarization_C = variables.Depolarization_C();
								double Depolarization_V = variables.Depolarization_V();
								double Depolarization_W = variables.Depolarization_W();

								// depolarization factors MC
								double mc_Depolarization_A = mc_variables.Depolarization_A();
								double mc_Depolarization_B = mc_variables.Depolarization_B();
								double mc_Depolarization_C = mc_variables.Depolarization_C();
								double mc_Depolarization_V = mc_variables.Depolarization_V();
								double mc_Depolarization_W = mc_variables.Depolarization_W();

								boolean matching_e = false;
								boolean matching_p1 = false;
								boolean matching_p2 = false;

								int matching_p1_pid = 0;
								int mc_p1_parent_index = 0;
								for (int current_part = 0; current_part < lundBank.rows(); current_part++) {
									int pid = lundBank.getInt("pid", current_part);
									if (matching_p1) { continue; }
									double mc_px = lundBank.getFloat("px", current_part);
									double mc_py = lundBank.getFloat("py", current_part);
									double mc_pz = lundBank.getFloat("pz", current_part);

									double mc_phi = phi_calculation(mc_px, mc_py);
									double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

									double exp_phi = phi_calculation(exp_p1.px(), exp_p1.py());
									double exp_theta = theta_calculation(exp_p1.px(), exp_p1.py(), 
										exp_p1.pz());

									matching_p1 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
										Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
									if (matching_p1) {
										matching_p1_pid = pid;
										mc_p1_parent_index = lundBank.getInt("parent", current_part)-1;
									}
								}
								int mc_p1_parent = lundBank.getInt("pid", mc_p1_parent_index);

								int matching_p2_pid = 0;
								int mc_p2_parent_index = 0;
								for (int current_part = 0; current_part < lundBank.rows(); current_part++) {
									int pid = lundBank.getInt("pid", current_part);
									if (matching_p2) { continue; }
									double mc_px = lundBank.getFloat("px", current_part);
									double mc_py = lundBank.getFloat("py", current_part);
									double mc_pz = lundBank.getFloat("pz", current_part);

									double mc_phi = phi_calculation(mc_px, mc_py);
									double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

									double exp_phi = phi_calculation(exp_p2.px(), exp_p2.py());
									double exp_theta = theta_calculation(exp_p2.px(), exp_p2.py(), 
										exp_p2.pz());

									matching_p2 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
										Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
									if (matching_p2) {
										matching_p2_pid = pid;
										mc_p2_parent_index = lundBank.getInt("parent", current_part)-1;
									} 
								}
								int mc_p2_parent = lundBank.getInt("pid", mc_p2_parent_index);

								matching_e = false;
								matching_p1 = false;
								matching_p2 = false;

								int matching_e_pid = 0;
								int mc_e_parent_index = 0;
								for (int current_part = 0; current_part < mcBank.rows(); current_part++) {
									int pid = mcBank.getInt("pid", current_part);
									if (matching_e) { continue; }
									double mc_px = mcBank.getFloat("px", current_part);
									double mc_py = mcBank.getFloat("py", current_part);
									double mc_pz = mcBank.getFloat("pz", current_part);

									double mc_phi = phi_calculation(mc_px, mc_py);
									double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

									double exp_phi = phi_calculation(exp_e.px(), exp_e.py());
									double exp_theta = theta_calculation(exp_e.px(), exp_e.py(), exp_e.pz());

									matching_e = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
										Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
									if (matching_e) {
										matching_e_pid = pid;
									}
								}

								matching_p1_pid = 0;
								mc_p1_parent_index = 0;
								for (int current_part = 0; current_part < mcBank.rows(); current_part++) {
									int pid = mcBank.getInt("pid", current_part);
									if (matching_p1) { continue; }
									double mc_px = mcBank.getFloat("px", current_part);
									double mc_py = mcBank.getFloat("py", current_part);
									double mc_pz = mcBank.getFloat("pz", current_part);

									double mc_phi = phi_calculation(mc_px, mc_py);
									double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

									double exp_phi = phi_calculation(exp_p1.px(), exp_p1.py());
									double exp_theta = theta_calculation(exp_p1.px(), exp_p1.py(), 
										exp_p1.pz());

									matching_p1 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
										Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
									if (matching_p1) {
										matching_p1_pid = pid;
									}
								}

								matching_p2_pid = 0;
								mc_p2_parent_index = 0;
								for (int current_part = 0; current_part < mcBank.rows(); current_part++) {
									int pid = mcBank.getInt("pid", current_part);
									if (matching_p2) { continue; }
									double mc_px = mcBank.getFloat("px", current_part);
									double mc_py = mcBank.getFloat("py", current_part);
									double mc_pz = mcBank.getFloat("pz", current_part);

									double mc_phi = phi_calculation(mc_px, mc_py);
									double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

									double exp_phi = phi_calculation(exp_p2.px(), exp_p2.py());
									double exp_theta = theta_calculation(exp_p2.px(), exp_p2.py(), 
										exp_p2.pz());

									matching_p2 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
										Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
									if (matching_p2) {
										matching_p2_pid = pid;
									} 
								}

								file.append(e_p+" "+mc_e_p+" "+e_theta+" "+mc_e_theta+" ");
								file.append(p1_p+" "+mc_p1_p+" "+p1_theta+" "+mc_p1_theta+" ");
								file.append(p2_p+" "+mc_p2_p+" "+p2_theta+" "+mc_p2_theta+" ");
								file.append(Q2+" "+mc_Q2+" "+W+" "+mc_W+" ");
								file.append(Mx+" "+mc_Mx+" "+Mx1+" "+mc_Mx1+" "+Mx2+" "+mc_Mx2+" ");
								file.append(x+" "+mc_x+" "+y+" "+mc_y+" "+z+" "+mc_z+" ");
								file.append(z1+" "+mc_z1+" "+z2+" "+mc_z2+" "+Mh+" "+mc_Mh+" ");
								file.append(xF+" "+mc_xF+" "+xF1+" "+mc_xF1+" "+xF2+" "+mc_xF2+" ");
								file.append(pT+" "+mc_pT+" "+pT1+" "+mc_pT1+" "+pT2+" "+mc_pT2+" ");
								file.append(pTpT+" "+mc_pTpT+" "+zeta+" "+mc_zeta+" ");
								file.append(eta+" "+mc_eta+" "+eta1+" "+mc_eta1+" "+eta2+" "+mc_eta2+" ");
								file.append(Delta_eta+" "+mc_Delta_eta+" ");
								file.append(eta1_gN+" "+mc_eta1_gN+" "+eta2_gN+" "+mc_eta2_gN+" ");
								file.append(phi1+" "+mc_phi1+" "+phi2+" "+mc_phi2+" ");
								file.append(Delta_phi+" "+mc_Delta_phi+" ");
								file.append(phih+" "+mc_phih+" "+phiR+" "+mc_phiR+" "+theta+" "+mc_theta+" ");
								file.append(matching_e_pid+" "+matching_p1_pid+" "+matching_p2_pid+" ");
								file.append(mc_p1_parent+" "+mc_p2_parent+"\n");

							}
						}
					}
				}
			}

			println(); println();
			print("1: e_p, 3: e_theta, 5: p1_p, 7: p1_theta, 9: p2_p, 11: p2_theta, ");
			print("13: Q2, 15: W, 17: Mx, 19: Mx1, 21: Mx2, 23: x, 25: y, 27: z ");
			print("29: z1, 31: z2, 33: Mh, 35: xF, 37: xF1, 39: xF2, ");
			print("41: pT, 43: pT1, 45: pT2, 47: pTpT, 49: zeta, ");
			print("51: eta, 53: eta1, 55: eta2, 57: Delta_eta, 59: eta1_gN, 61: eta2_gN, ");
			print("63: phi1, 65: phi2, 67: Delta_phi, 69: phih, 71: phiR, 73: theta ");
			print("74: matching e id, 75: matching p1 pid, 76: matching p2 pid ");
			print("77: p1 parent pid, 78: p2 parent pid\n");

			println(); println();
			println("Set p1 PID = "+p1_Str+"\n");
			println("Set p2 PID = "+p2_Str+"\n");
			println("output file is: "+file);
			println();
		}
	}
}