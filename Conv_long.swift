import "nek";
type file;

/* Inputs: standard dictionary, user file, and list of viscosities */
string prefix="cnv_long";
file json <"Conv_long.json">;
file tusr <"single_mode.tusr">;

int max_dof = 128;
//int max_dof = 256;
//int max_dof = 512;
//float[] orders = [4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0];
//int[] orders = [4, 6, 8, 10, 12, 14, 16, 32];
//int[] orders = [4, 8, 16, 32];
int[] orders = [8];
//float[] elms = [64.0, 40.0, 32.0, 26.0, 20.0, 18.0, 16.0];
//int[] elms = [4, 8, 12, 16, 20, 24, 28, 32];
//int[] elms = [4, 8, 16, 32, 64];
int[] elms = [16];

int[int] dealias;
dealias[4]  =  6;
dealias[6]  = 10;
dealias[8]  = 12;
dealias[10] = 16;
dealias[12] = 18;
dealias[14] = 22;
dealias[16] = 24;
dealias[32] = 48;

int job_wall = 60;

int nstep = 24;
int job_step = 1;
int io_step = 1;

float io_time = 1.0;
float job_time = 4.0;
int outputs_per_job = 4;

int j0 = 0;
string analysis = "RTI";
int post_nodes = 0;

int aspect = 32;
//int aspect = 4;

string cwd = arg("cwd", ".");

foreach order,i in orders {
  foreach elm,ii in elms {
    if (order * elm <= max_dof) {

      int mode;
      if (order <= 16){
        mode = 64;
      } else {
        mode = 32;
      }

      int ncore;
      if (toInt(order*order*order*elm*elm*elm*aspect /16384) < 1){
        ncore = 1;
      } else if (toInt(order*order*order*elm*elm*elm*aspect /16384) > elm*elm*elm*aspect){
        ncore = elm*elm*elm*aspect;
      } else {
        ncore = toInt(order*order*order*elm*elm*elm*aspect / 16384);
      }

      int nodes = toInt(ncore / mode);
      int nwrite;
      if (nodes < 4) {
        nwrite = 1;
      }else{
        nwrite = toInt(nodes / 4);
      }

     float dt_max = (2/(elm*(order-1)*(order-1)))/0.0558519;
 
      // Pick a directory to run in 
      string tdir = sprintf("./%s_o%i_e%i", prefix, order, elm);
      string tdir_f = sprintf("%s/%s_o%i_e%i", cwd, prefix, order, elm);
      string name = sprintf("./%s_o%i_e%i", prefix, order, elm);
      file foo <single_file_mapper; file=strcat("mkdirs/", tdir, ".out")>;
      (foo) = mkdir(tdir_f);
  
      // Construct input files and build the nek5000 executable 
      file base     <single_file_mapper; file=sprintf("%s/%s.json", tdir, name)>;
      file rea      <single_file_mapper; file=sprintf("%s/%s.rea",  tdir, name)>;
      file map      <single_file_mapper; file=sprintf("%s/%s.map",  tdir, name)>;
      file usr      <single_file_mapper; file=sprintf("%s/%s.usr",  tdir, name)>;
      file size_mod <single_file_mapper; file=sprintf("%s/size_mod.F90",  tdir)>;
 
      string override = sprintf("{ \"order\": %i, \"dealiasing_order\": %i, \"shape_mesh\": [%i, %i, %i], \"procs\": %i, \"io_files\": %i, \"dt\": %f, \"io_time\": %f }", 
                                order, dealias[order], elm, elm, elm*aspect, ncore, -nwrite, dt_max, io_time);
       
      (usr, rea, map, base, size_mod) = genrun_str (json, tusr, name, tdir_f, override, foo);
      file nek5000 <single_file_mapper; file=sprintf("%s/nek5000", tdir_f, name)>;
      (nek5000) = makenek(tdir_f, "/projects/HighAspectRTI/nek/", name, usr, size_mod);

      int ret_l;
      (ret_l) =  series(nek5000, nodes, mode,
                        cwd, tdir, tdir_f, 
                        name, nwrite, base, usr, job_time, nstep, job_step, outputs_per_job, 
                        analysis=analysis, job_wall=job_wall, post_nodes=post_nodes, j0=j0);
    }
  }
}
