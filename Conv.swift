import "nek";
type file;

/* Inputs: standard dictionary, user file, and list of viscosities */
string prefix="cnv";
file json <"Conv.json">;
file tusr <"single_mode.tusr">;

string pname="order";
string qname="courant";
//float[] pvals = [4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0];
int[] pvals = [4, 8];
//float[] qvals = [64.0, 40.0, 32.0, 26.0, 20.0, 18.0, 16.0];
int[] qvals = [32, 16];

int[int] dealias;
dealias[4]  =  6;
dealias[6]  = 10;
dealias[8]  = 12;
dealias[10] = 16;
dealias[12] = 18;
dealias[14] = 22;
dealias[16] = 24;

int nwrite=128;

int nodes = 128;
int mode = 64;
int job_wall = 60;

int nstep = 48;
int io_step = 1;
int job_step = 1;

float io_time = 0.5;
float job_time = 1.0;

int j0 = 0;
string analysis = "RTI";
boolean legacy = false;
int post_nodes = 0;


int aspect = 8;

int outputs_per_job;
if (job_time > 0.0 && io_time > 0.0){
  trace("Foo");
  //outputs_per_job = toInt(job_time / io_time);
  outputs_per_job = 2; //toInt(sprintf("%f", job_time / io_time));
  //tracef("Foo!");
  //outputs_per_job = 2;
} else {
  outputs_per_job = job_step %/ io_step;
}
string cwd = arg("cwd", ".");

  foreach pval,i in pvals {
    foreach qval,ii in qvals {
      if (pval * qval <= 256) {
 
        /* Pick a directory to run in */
        string tdir = sprintf("./%s_o%i_e%i", prefix, pval, qval);
        string name = sprintf("./%s_o%i_e%i", prefix, pval, qval);
        file tdir_f  <single_file_mapper; file=strcat(cwd,"/",tdir)>;
        file foo     <single_file_mapper; file=strcat(cwd,"/",tdir, "/foo")>;
        (foo) = mkdir(tdir_f);
    
        /* Construct input files and build the nek5000 executable */
        file base     <single_file_mapper; file=sprintf("%s/%s.json", tdir, name)>;
        file rea      <single_file_mapper; file=sprintf("%s/%s.rea",  tdir, name)>;
        file map      <single_file_mapper; file=sprintf("%s/%s.map",  tdir, name)>;
        file usr      <single_file_mapper; file=sprintf("%s/%s.usr",  tdir, name)>;
        //file size_mod <single_file_mapper; file=sprintf("%s/SIZE",  tdir, name)>;
        file size_mod <single_file_mapper; file=sprintf("%s/size_mod.F90",  tdir)>;
   
        string override = sprintf("{ \"order\": %i, \"dealiasing_order\": %i, \"shape_mesh\": [%i, %i, %i]}", pval, dealias[pval], qval, qval, qval*aspect);
         
        (usr, rea, map, base, size_mod) = genrun_str (json, tusr, name, tdir_f, override, foo, _legacy=legacy);
        file nek5000 <single_file_mapper; file=sprintf("%s/nek5000", tdir, name)>;
        (nek5000) = makenek(tdir_f, "/projects/HighAspectRTI/nek/", name, usr, size_mod, _legacy=legacy);
  
        int ret_l;
        (ret_l) =  series(nek5000, nodes, mode,
                          cwd, tdir, tdir_f, 
                          name, nwrite, base, usr, job_time, nstep, job_step, outputs_per_job, 
                          analysis=analysis, job_wall=job_wall, post_nodes=post_nodes, j0=j0);
      }
    }
  }

