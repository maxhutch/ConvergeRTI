import "nek";
type file;

/* Inputs: standard dictionary, user file, and list of viscosities */
type config {
  int max_dof;
  int min_dof;
  int[] orders;
  int[] elms;
  int nstep;
  int job_step;
  int io_step;
  float io_time;
  float job_time;
  int outputs_per_job;
  int aspect;
  int j0;
  string prefix;
  string json_name;
  string tusr_name;
}

config c;
c = readStructured(arg("in"));
string prefix    = strcut(c.prefix, "\"([^ ]*)\"");
string json_name = strcut(c.json_name, "\"([^ ]*)\"");
string tusr_name = strcut(c.tusr_name, "\"([^ ]*)\"");
file json <single_file_mapper;file=json_name>;
file tusr <single_file_mapper;file=tusr_name>;


int[int] dealias;
dealias[4]  =  6;
dealias[6]  = 10;
dealias[8]  = 12;
dealias[10] = 16;
dealias[12] = 18;
dealias[14] = 22;
dealias[16] = 24;
dealias[32] = 48;


int[int] nelm_per_proc;
nelm_per_proc[32] = 1; // note that mode is changing for this case
nelm_per_proc[16] = 4;
nelm_per_proc[14] = 4;
nelm_per_proc[12] = 8;
nelm_per_proc[10] = 16;
nelm_per_proc[8]  = 32;
nelm_per_proc[6]  = 64;
nelm_per_proc[4]  = 256;

int job_wall = 60;

string analysis = "RTI";
int post_nodes = 0;


string cwd = arg("cwd", ".");

foreach order,i in c.orders {
  foreach elm,ii in c.elms {
    if (order * elm <= c.max_dof && order * elm >= c.min_dof) {

      int mode;
      if (order <= 16){
        mode = 64;
      } else {
        mode = 32;
      }

      int ncore;
      if (toInt(elm*elm*elm*c.aspect / nelm_per_proc[order]) < 1){
        ncore = 1;
      } else if (toInt(elm*elm*elm*c.aspect /nelm_per_proc[order]) > elm*elm*elm*c.aspect){
        ncore = elm*elm*elm*c.aspect;
      } else {
        ncore = toInt(elm*elm*elm*c.aspect / nelm_per_proc[order]);
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
 
      string override = sprintf(
                        "{ \"order\": %i, \"dealiasing_order\": %i, \"shape_mesh\": [%i, %i, %i], " 
                      + "\"procs\": %i, \"io_files\": %i, \"dt\": %f, \"io_time\": %f, "
                      + "\"root_mesh\": [%f, %f, %f], \"extent_mesh\": [%f, %f, %f]}", 
                                order, dealias[order], elm, elm, elm*c.aspect, 
                                ncore, -nwrite, dt_max, c.io_time, 
                                0.0, 0.0, -c.aspect/4, 0.5, 0.5, c.aspect/4);
       
      (usr, rea, map, base, size_mod) = genrun_str (json, tusr, name, tdir_f, override, foo);
      file nek5000 <single_file_mapper; file=sprintf("%s/nek5000", tdir_f, name)>;
      (nek5000) = makenek(tdir_f, "/projects/HighAspectRTI/nek/", name, usr, size_mod);

      int ret_l;
      (ret_l) =  series(nek5000, nodes, mode,
                        cwd, tdir, tdir_f, 
                        name, nwrite, base, usr, c.job_time, c.nstep, c.job_step, c.outputs_per_job, 
                        analysis=analysis, job_wall=job_wall, post_nodes=post_nodes, j0=c.j0);
    }
  }
}
