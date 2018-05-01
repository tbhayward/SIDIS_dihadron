/**
 *
 * @author Timothy B. Hayward
 */

package SIDIS_fitter;

import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

public class SIDIS_fitter extends GenericKinematicFitter {

    protected final Double mybeam;

    public SIDIS_fitter(double beam) {
        super(beam);
        mybeam = beam;
    } 
    /**
     * Returns PhysicsEvent object with reconstructed particles.
     *
     * @param event - DataEvent object
     * @return PhysicsEvent : event containing particles.
     */
    @Override
    public PhysicsEvent getPhysicsEvent(DataEvent event) {
        
        boolean banks_test = true; // check to see if the event has all of the banks present
        if (!(event.hasBank("REC::Particle"))) {
            banks_test = false;
        } else if (!(event.hasBank("REC::Cherenkov"))) {
            banks_test = false;
        }
        
        if (banks_test) {
            PhysicsEvent physEvent = new PhysicsEvent();
            HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Particle");
            HipoDataBank ccpbBank = (HipoDataBank) event.getBank("REC::Cherenkov"); // load cherenkov counter bank
            
            for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
                float vx = eventBank.getFloat("vx",current_part);
                float vy = eventBank.getFloat("vy",current_part);
                float vz = eventBank.getFloat("vz",current_part);
                boolean vertex_test = true; // check to see if vertex within 5cm of origin
                float vertex_cutoff = 4;
                if (Math.abs(vx)>vertex_cutoff) {
                    vertex_test = false;
                } else if (Math.abs(vy)>vertex_cutoff) {
                    vertex_test = false;
                } else if (Math.abs(vz)>vertex_cutoff) {
                    vertex_test = false;
                }
                if (vertex_test) {
                    
                    short nphe = -1; // number of photoelectrons in CC, set to -1 to allow for "0" as a test for hadrons
                    for (int cc_hit = 0; cc_hit < ccpbBank.rows(); cc_hit++) {
                        int cc_pindex = ccpbBank.getShort("pindex", cc_hit); // row # in particle bank associated
                        if (cc_pindex==current_part) { // if the current 
                            nphe = ccpbBank.getShort("nphe", cc_hit);
                        }
                    }
                    float px = eventBank.getFloat("px", current_part);
                    float py = eventBank.getFloat("py", current_part);
                    float pz = eventBank.getFloat("pz", current_part);
                    if (nphe > 2) { // ultrarelativistic leptons produce larger signal in CC than hadrons
                        // adds lepton, negative the charge because PID of electron = 11
                        Particle part = new Particle(-11*eventBank.getByte("charge", current_part),px,py,pz,vx,vy,vz);
                        physEvent.addParticle(part);
                    } else if (nphe > -1) {
                        // adds pion
                        Particle part = new Particle(211*eventBank.getByte("charge", current_part),px,py,pz,vx,vy,vz);
                        physEvent.addParticle(part);
                    }
                }
            }
          
            return physEvent;
        }
    return new PhysicsEvent(this.mybeam);
    }
}
    
