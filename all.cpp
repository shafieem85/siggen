#include <pybind11/pybind11.h>
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
// #include <unistd.h>
#include <io.h>  
// <io.h>  similar to <unistd.h> but for windows
#include <string.h>
#include <math.h>
#include <time.h>
// #include "mjd_siggen.h"
// #include "detector_geometry.h"
#include "all.h"


/* ************************* point.h *************************** */
#define SQ(x) ((x)*(x))

/* norm
 * returns the absolute value of vector v
 * e is unit vector in v's direction
 */
float vector_norm(vector v, vector *e){
  float length;
  
  length = vector_length(v);

  e->x = v.x / length;
  e->y = v.y / length;
  e->z = v.z / length;
  
  return length;
}

/* vector_length
 * returns the length of vector v
 */
float vector_length(vector v){
  return sqrt(SQ(v.x) + SQ(v.y) + SQ(v.z));
}

/* distance
 * returns distance between points
 */
float distance(point pt1, point pt2){
  float d;

  d = sqrt(SQ(pt1.x - pt2.x) + SQ(pt1.y - pt2.y) + SQ(pt1.z - pt2.z));

  return d;
}

/* vector_add */
vector vector_add(vector v1, vector v2){
  vector v;
  
  v.x = v1.x + v2.x;
  v.y = v1.y + v2.y;
  v.z = v1.z + v2.z;

  return v;
}

/*vector subtraction, v1-v2*/
vector vector_sub(vector v1, vector v2){
  vector v;
  
  v.x = v1.x - v2.x;
  v.y = v1.y - v2.y;
  v.z = v1.z - v2.z;
  
  return v;
}

/*vector dot product, v1*v1*/
float dot_prod(vector v1, vector v2){
  return  v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

/*vector cross product, v1 X v2*/
vector cross_prod(vector v1, vector v2){
  vector v;
  
  v.x = v1.y*v2.z - v1.z*v2.y;
  v.y = v1.z*v2.x - v1.x*v2.z;
  v.z = v1.x*v2.y - v1.y*v2.x;

  return v;
}

/*scale vector by factor -- better name?*/
vector vector_scale(vector v, float factor){
  vector res;

  res.x = v.x*factor;
  res.y = v.y*factor;
  res.z = v.z*factor;

  return res;
}

/*rotate vector by angle_deg degrees*/
vector vector_rotate_z(vector v, float angle_deg){
  float angle_rad;
  vector res;
  angle_rad = angle_deg/180.0*M_PI;

  res.x = cos(angle_rad)*v.x - sin(angle_rad)*v.y;
  res.y = sin(angle_rad)*v.x + cos(angle_rad)*v.y;
  res.z = v.z;

  return res;
}

/* pt_to_str
 * string representation of point pt
 * written to string str, which has guaranteed length len 
 */
char *pt_to_str(char *str, int len, point pt){
  snprintf(str, len, "(%.1f %.1f %.1f)", pt.x, pt.y, pt.z);
  return str;
}

#undef SQ

/* ************************* cyl_point.h *************************** */
struct cyl_pt cart_to_cyl(struct point cart_pt){
  struct cyl_pt cyl_pt;

  cyl_pt.r = sqrt(cart_pt.x*cart_pt.x + cart_pt.y*cart_pt.y);
  if (cart_pt.x == 0.0){
    if (cart_pt.y > 0.0) cyl_pt.phi = M_PI/2;
    else cyl_pt.phi = -M_PI/2;
  }else{
    cyl_pt.phi = atan(cart_pt.y/cart_pt.x);
    if (cart_pt.x < 0) cyl_pt.phi += M_PI;
  }
  cyl_pt.z = cart_pt.z;
  
  return cyl_pt;
}

struct point cyl_to_cart(struct cyl_pt cyl_pt){
  struct point cart_pt;

  cart_pt.x = cyl_pt.r * cos(cyl_pt.phi);
  cart_pt.y = cyl_pt.r * sin(cyl_pt.phi);
  cart_pt.z = cyl_pt.z;

  return cart_pt;
}

#define SQUARED(x) ((x)*(x))
float cart_distance(struct point pt1, struct point pt2){
  return sqrt(SQUARED(pt1.x-pt2.x) + SQUARED(pt1.y-pt2.y) 
	      + SQUARED(pt1.z-pt2.z));
}
#undef SQUARED


float cyl_distance(struct cyl_pt pt1, struct cyl_pt pt2){
  /*yes, I'm very lazy. It's one of my best traits*/
  return cart_distance(cyl_to_cart(pt1), cyl_to_cart(pt2));
}

struct cyl_pt cyl_diff(struct cyl_pt from, struct cyl_pt to){
  struct cyl_pt v;

  v.r = from.r - to.r;
  v.phi = from.phi - to.phi;
  v.z = from.z - to.z;

  return v;
}

float vector_norm_cyl(struct cyl_pt pt, struct cyl_pt *e){
  float norm;

  norm = sqrt(pt.r*pt.r + pt.z*pt.z);
  e->phi = pt.phi;
  e->r = pt.r/norm;
  e->z = pt.z/norm;

  return norm;
}


#define SQ(x) ((x)*(x))
/* outside_detector
   returns 1 if pt is outside the detector, 0 if inside detector
*/
int outside_detector(point pt, MJD_Siggen_Setup *setup){
  float r, z, r1, z1, br, a, b;

  z = pt.z;
  if (z > setup->zmax || z < 0) return 1;

  r = sqrt(SQ(pt.x)+SQ(pt.y));
  if (r > setup->rmax) return 1;
  r1 = setup->rmax - r;  // distance from outer radius
  z1 = setup->zmax - z;  // distance from top of crystal

  /* check point contact */
  if (z < setup->pc_length && r < setup->pc_radius) {
    if (!setup->bulletize_PC) return 1;
    if (setup->pc_length > setup->pc_radius) {
      a = setup->pc_length - setup->pc_radius;
      if (z < a || SQ(z-a) + SQ(r) < SQ(setup->pc_radius)) return 1;
    } else {
      a = setup->pc_radius - setup->pc_length;
      if (r < a || SQ(z) + SQ(r-a) < SQ(setup->pc_length)) return 1;
    }
    return 0;
  }

  /* check ditch */
  if (z < setup->ditch_depth  &&
      setup->ditch_thickness > 0 && setup->wrap_around_radius > 0 &&
      r < setup->wrap_around_radius &&
      r > setup->wrap_around_radius - setup->ditch_thickness) return 1;

  /* check hole */
  if ( r < setup->hole_radius &&
      z1 < setup->hole_length) {
    b = setup->zmax - setup->hole_length + setup->hole_bullet_radius;
    if (z > b) return 1;
    a = setup->hole_radius - setup->hole_bullet_radius;
    if (r < a || SQ(b-z) + SQ(r-a) < SQ(setup->hole_bullet_radius)) return 1;
  }

  /* check inner taper of hole */
  if (z1 < setup->inner_taper_length &&
      r  < setup->hole_radius +
            ((setup->inner_taper_length - z1) *
              setup->inner_taper_width / setup->inner_taper_length)) return 1;      
  /* check outer taper of crystal */
  if (z1 < setup->outer_taper_length &&
      r1 < ((setup->outer_taper_length - z1) *
              setup->outer_taper_width / setup->outer_taper_length)) return 1;
  /* check 45-degree bottom outer taper of crystal */
  if ( z < setup->bottom_taper_length &&
      r1 < setup->bottom_taper_length - z) return 1;

  /* check bulletizations */
  br = setup->top_bullet_radius;
  // adjust top bulletization position for top outer taper
  a = 0;
  if (setup->outer_taper_length > br)
    a = ((setup->outer_taper_length - br) *
         setup->outer_taper_width / setup->outer_taper_length);
  if (z1 < br &&
      r1 < br + a &&
      SQ(br - r1 + a) + SQ(br - z1) > br*br) return 1;
  br = setup->bottom_bullet_radius;
  if ( z < br &&
      r1 < br &&
      SQ(br - r1) + SQ(br - z ) > br*br) return 1;

  return 0;
}

/* ************************* detector_geometry.h *************************** */
int outside_detector_cyl(cyl_pt pt, MJD_Siggen_Setup *setup){
  float r, z, r1, z1, br, a, b;

  z = pt.z;
  if (z > setup->zmax || z < 0) return 1;

  r = pt.r;
  if (r < 0) r = -r;
  if (r > setup->rmax) return 1;
  r1 = setup->rmax - r;  // distance from outer radius
  z1 = setup->zmax - z;  // distance from top of crystal

  /* check point contact */
  if (z < setup->pc_length && r < setup->pc_radius) {
    if (!setup->bulletize_PC) return 1;
    if (setup->pc_length > setup->pc_radius) {
      a = setup->pc_length - setup->pc_radius;
      if (z < a || SQ(z-a) + SQ(r) < SQ(setup->pc_radius)) return 1;
    } else {
      a = setup->pc_radius - setup->pc_length;
      if (r < a || SQ(z) + SQ(r-a) < SQ(setup->pc_length)) return 1;
    }
    return 0;
  }
  /* check ditch */
  if (z <= setup->ditch_depth  &&
      setup->ditch_thickness > 0 && setup->wrap_around_radius > 0 &&
      r < setup->wrap_around_radius &&
      r > setup->wrap_around_radius - setup->ditch_thickness) return 1;

  /* check hole */
  if ( r < setup->hole_radius &&
      z1 < setup->hole_length) {
    b = setup->zmax - setup->hole_length + setup->hole_bullet_radius;
    if (z > b) return 1;
    a = setup->hole_radius - setup->hole_bullet_radius;
    if (r < a || SQ(b-z) + SQ(r-a) < SQ(setup->hole_bullet_radius)) return 1;
  }
  /* check inner taper of hole */
  if (z1 < setup->inner_taper_length &&
      r  < setup->hole_radius +
            ((setup->inner_taper_length - z1) *
              setup->inner_taper_width / setup->inner_taper_length)) return 1;      
  /* check outer taper of crystal */
  if (z1 < setup->outer_taper_length &&
      r1 < ((setup->outer_taper_length - z1) *
              setup->outer_taper_width / setup->outer_taper_length)) return 1;
  /* check 45-degree bottom outer taper of crystal */
  if ( z < setup->bottom_taper_length &&
      r1 < setup->bottom_taper_length - z) return 1;

  /* check bulletizations */
  br = setup->top_bullet_radius;
  // adjust top bulletization position for top outer taper
  a = 0;
  if (setup->outer_taper_length > br)
    a = ((setup->outer_taper_length - br) *
         setup->outer_taper_width / setup->outer_taper_length);
  if (z1 < br &&
      r1 < br + a &&
      SQ(br - r1 + a) + SQ(br - z1) > br*br) return 1;
  br = setup->bottom_bullet_radius;
  a = setup->Li_thickness;   // FIXME ? added for fieldgen
  if ( z < br + a &&
      r1 < br &&
      SQ(br - r1) + SQ(br - z + a) > br*br) return 1;

  return 0;
}
#undef SQ

/* ************************* read_config.h *************************** */

int read_config(char *config_file_name, MJD_Siggen_Setup *setup) {

  /* reads and parses configuration file of name config_file_name
     fills in values of MJD_Siggen_Setup setup, defined in mjd_siggen.h
     returns 0 on success, 1 otherwise
  */

  char key_word[][32] = {
    "xtal_length",
    "xtal_radius",
    "top_bullet_radius",
    "bottom_bullet_radius",
    "pc_length",
    "pc_radius",
    "bulletize_PC",
    "bottom_taper_length",  // note: these two keywords must stay adjacent, in this order
    "taper_length",         // <- for compatibility with old mjd config files, == "bottom taper_length
    "outer_taper_length",
    "outer_taper_width",
    "taper_angle",
    "inner_taper_length",
    "inner_taper_width",
    "hole_length_gap",     // can use gap (i.e. xtal_length - hole_length) instead. This keyword must be before "hole_length".
    "hole_length",
    "hole_radius",
    "hole_bullet_radius",
    "wrap_around_radius",
    "ditch_depth",
    "ditch_thickness",
    "Li_thickness",
    "vacuum_gap",
    "xtal_grid",
    "impurity_z0",
    "impurity_gradient",
    "impurity_quadratic",
    "impurity_surface",
    "impurity_radial_add",
    "impurity_radial_mult",
    "impurity_rpower",
    "xtal_HV",
    "surface_drift_vel_factor",

    "drift_name",
    "field_name",
    "wp_name",
    "xtal_temp",
    "preamp_tau",
    "time_steps_calc",
    "step_time_calc",
    "step_time_out",
    "charge_cloud_size",
    "use_diffusion",
    "energy",
    "verbosity_level",
    "max_iterations",
    "write_field",
    "write_WP",
    ""
  };

  int   ii, i, l, n=0, ok, iint = 0, diameter = 0;
  float fi;
  char  *c, line[256], line2[256], name[256];
  FILE  *file;


  /* initialize everything to zero... */
  memset(setup, 0, sizeof(*setup));
  /* ...except for impurity_radial_mult */
  setup->impurity_radial_mult = 1.0f;      // 1.0 is neutral (no radial gradient)
  setup->surface_drift_vel_factor = 1.0f;  // normal fast drift, same as bulk

  if (!(file = fopen(config_file_name, "r"))) {
    printf("\nERROR: config file %s does not exist?\n", config_file_name);
    return 1;
  }
  /* read config file */
  printf("\nReading values from config file %s\n", config_file_name);
  while (fgets(line, sizeof(line), file)) {
    n++;
    /* ignore comments and blank lines */
    if (strlen(line) < 3 || *line == ' ' || *line == '\t' || *line == '#') continue;
    /* if line contains "_diam" (for diameter) replace with "_radius"
       this allows the user to specify diametwrs instead of radii */
    diameter = 0;
    if ((c = strstr(line, "_diam"))) {
      diameter = 1;
      if (setup->verbosity >= CHATTY) printf("Line: %s", line);
      strcpy(line2, c+5);
      strcpy(c, "_radius");
      strcpy(c+7, line2);
      if (setup->verbosity >= CHATTY) printf("   -->  %s", line);      
    }

    for (i=0; (l=strlen(key_word[i])) > 0; i++) {
      if (!strncmp(line, key_word[i], l)) {
	/* line starts with key_word[i] */
      printf("Line: %s", line);
	if (line[l] != ' ' && line[l] != '\t') {
	  ok = 0;
	} else {
	  // for (c = line + l; *c != ' ' && *c != '\t'; c++) ;
	  /* find next non-white-space char */
	  for (c = line + l; *c == ' ' || *c == '\t'; c++) ;
	  name[0] = 0;
	  ii = iint = 0;
	  fi = 0;
	  if (strstr(key_word[i], "_name")) {
	    /* extract character string for file name */
	    for (ok=0; ok<256 && *c != ' ' && *c != '\t' && *c != '\n' &&  *c != '\r'; ok++) {
	      name[ok] = *c;
	      c++;
	    }
	    name[ok] = '\0';	    
	  } else if (!strncmp("time_steps_calc", key_word[i], l) ||
		     !strncmp("use_diffusion", key_word[i], l) ||
		     !strncmp("verbosity_level", key_word[i], l) ||
		     !strncmp("max_iterations", key_word[i], l) ||
		     !strncmp("write_field", key_word[i], l) ||
		     !strncmp("write_WP", key_word[i], l) ||
		     !strncmp("bulletize_PC", key_word[i], l)) {
	    /* extract integer value */
	    ok = sscanf(c, "%d", &ii);
	    iint = 1;
	  } else {
	    /* extract float value */
	    ok = sscanf(c, "%f", &fi);
            if (diameter) fi /= 2.0;
	  }
	}
	if (ok < 1) {
	  printf("ERROR reading %s from config file %s\n"
		 "   ...line number %d is: %s",
		 key_word[i], config_file_name, n, line);
	  return 1;
	}

        if (!strncmp(key_word[i], "taper_length", l)) {
          i--;                      // use previous keyword = "bottom_taper_length" instead
          l=strlen(key_word[i]);    // for backwards compatibiity wth old mjd config files
        }

	if (strstr(key_word[i], "verbosity_level")) {
	  setup->verbosity = ii;
	} else if (strstr(key_word[i], "xtal_length")) {
	  setup->xtal_length = fi;
	} else if (strstr(key_word[i], "xtal_radius")) {
	  setup->xtal_radius = fi;
	} else if (strstr(key_word[i], "top_bullet_radius")) {
	  setup->top_bullet_radius = fi;
	} else if (strstr(key_word[i], "bottom_bullet_radius")) {
	  setup->bottom_bullet_radius = fi;
	} else if (strstr(key_word[i], "pc_length")) {
	  setup->pc_length = fi;
	} else if (strstr(key_word[i], "pc_radius")) {
	  setup->pc_radius = fi;
	} else if (strstr(key_word[i], "bulletize_PC")) {
	  setup->bulletize_PC = ii;
	} else if (strstr(key_word[i], "bottom_taper_length")) {
	  setup->bottom_taper_length = fi;
	} else if (strstr(key_word[i], "outer_taper_length")) {
	  setup->outer_taper_length = fi;
	} else if (strstr(key_word[i], "outer_taper_width")) {
	  setup->outer_taper_width = fi;
	} else if (strstr(key_word[i], "taper_angle")) {
	  setup->taper_angle = fi;
	} else if (strstr(key_word[i], "inner_taper_length")) {
	  setup->inner_taper_length = fi;
	} else if (strstr(key_word[i], "inner_taper_width")) {
	  setup->inner_taper_width = fi;
	} else if (strstr(key_word[i], "hole_length")) {
	  setup->hole_length = fi;
          /* the user can specify "hole_length_gap" = xtal_length - hole_length, instead of "hole_length" */
          if (strstr(line, "hole_length_gap")  && fi < setup->xtal_length) {
            setup->hole_length = setup->xtal_length - fi;
            if (setup->verbosity >= CHATTY) printf("   -->  hole_length: %f\n", setup->hole_length);
          }
	} else if (strstr(key_word[i], "hole_radius")) {
	  setup->hole_radius = fi;
	} else if (strstr(key_word[i], "hole_bullet_radius")) {
	  setup->hole_bullet_radius = fi;
	} else if (strstr(key_word[i], "wrap_around_radius")) {
	  setup->wrap_around_radius = fi;
	} else if (strstr(key_word[i], "ditch_depth")) {
	  setup->ditch_depth = fi;
	} else if (strstr(key_word[i], "ditch_thickness")) {
	  setup->ditch_thickness = fi;
	} else if (strstr(key_word[i], "Li_thickness")) {
	  setup->Li_thickness = fi;
	} else if (strstr(key_word[i], "vacuum_gap")) {
	  setup->vacuum_gap = fi;
	} else if (strstr(key_word[i], "xtal_grid")) {
	  setup->xtal_grid = fi;
	} else if (strstr(key_word[i], "impurity_z0")) {
	  setup->impurity_z0 = fi;
	} else if (strstr(key_word[i], "impurity_gradient")) {
	  setup->impurity_gradient = fi;
	} else if (strstr(key_word[i], "impurity_quadratic")) {
	  setup->impurity_quadratic = fi;
	} else if (strstr(key_word[i], "impurity_surface")) {
	  setup->impurity_surface = fi;
	} else if (strstr(key_word[i], "impurity_radial_add")) {
	  setup->impurity_radial_add = fi;
	} else if (strstr(key_word[i], "impurity_radial_mult")) {
	  setup->impurity_radial_mult = fi;
	} else if (strstr(key_word[i], "impurity_rpower")) {
	  setup->impurity_rpower = fi;
	} else if (strstr(key_word[i], "xtal_HV")) {
	  setup->xtal_HV = fi;
	} else if (strstr(key_word[i], "surface_drift_vel_factor")) {
	  setup->surface_drift_vel_factor = fi;
	} else if (strstr(key_word[i], "drift_name")) {
	  strncpy(setup->drift_name, name, 256);
	} else if (strstr(key_word[i], "field_name")) {
	  strncpy(setup->field_name, name, 256);
	} else if (strstr(key_word[i], "wp_name")) {
	  strncpy(setup->wp_name, name, 256);
	} else if (strstr(key_word[i], "xtal_temp")) {
	  setup->xtal_temp = fi;
	} else if (strstr(key_word[i], "preamp_tau")) {
	  setup->preamp_tau = fi;
	} else if (strstr(key_word[i], "time_steps_calc")) {
	  setup->time_steps_calc = ii;
	} else if (strstr(key_word[i], "step_time_calc")) {
	  setup->step_time_calc = fi;
	} else if (strstr(key_word[i], "step_time_out")) {
	  setup->step_time_out = fi;
	} else if (strstr(key_word[i], "charge_cloud_size")) {
	  setup->charge_cloud_size = fi;
	} else if (strstr(key_word[i], "use_diffusion")) {
	  setup->use_diffusion = ii;
	} else if (strstr(key_word[i], "energy")) {
	  setup->energy = fi;
	} else if (strstr(key_word[i], "max_iterations")) {
	  setup->max_iterations = ii;
	} else if (strstr(key_word[i], "write_field")) {
	  setup->write_field = ii;
	} else if (strstr(key_word[i], "write_WP")) {
	  setup->write_WP = ii;
	} else {
	  printf("ERROR; unrecognized keyword %s\n", key_word[i]);
	  return 1;
	}

	if (setup->verbosity >= CHATTY) {
	  // printf("%s", line);
	  if (iint) {
	    printf("%s: %d\n", key_word[i], ii);
	  } else if (strlen(name) > 0) {
	    printf("%s: %s\n", key_word[i], name);
	  } else {
	    printf("%s: %f\n", key_word[i], fi);
	  }
	}
	break;
      }
    }
  }
  fclose(file);

  if (setup->taper_angle > 0) {
    /* convert taper angle to taper widths */
    if (setup->outer_taper_length > 0) {
      setup->outer_taper_width =
        setup->outer_taper_length * tan(setup->taper_angle * 3.14159/180.0);
      printf("  ->>  outer taper width: %f\n", setup->outer_taper_width);
    }
    if (setup->inner_taper_length > 0) {
      setup->inner_taper_width =
        setup->inner_taper_length * tan(setup->taper_angle * 3.14159/180.0);
      printf("  ->>  inner taper width: %f\n", setup->inner_taper_width);
    }
  } else {
    /* convert taper width to taper angle */
    if (setup->outer_taper_length > 0 &&
        setup->outer_taper_width > 0)
      setup->taper_angle =
        atan(setup->outer_taper_width/setup->outer_taper_length) * 180.0/3.14159;
    else if (setup->inner_taper_length > 0 &&
             setup->inner_taper_width > 0)
      setup->taper_angle =
        atan(setup->inner_taper_width/setup->inner_taper_length) * 180.0/3.14159;
    if (setup->taper_angle > 0)
      printf("  ->>  taper angle: %f\n", setup->taper_angle);
  }
  if (setup->wrap_around_radius == 0 ||
      setup->wrap_around_radius > setup->xtal_radius - setup->bottom_taper_length)
    setup->wrap_around_radius = setup->xtal_radius - setup->bottom_taper_length;

  /* some consistency checks */
  if (setup->inner_taper_length > setup->hole_length - setup->hole_bullet_radius)
    setup->inner_taper_length = setup->hole_length - setup->hole_bullet_radius;
  if (setup->hole_bullet_radius > setup->hole_radius)
    setup->hole_bullet_radius = setup->hole_radius;
  if (setup->outer_taper_length > setup->xtal_length ||
      setup->inner_taper_length > setup->hole_length ||
      setup->hole_length > setup->xtal_length ||
      setup->hole_length < setup->hole_bullet_radius ||
      setup->inner_taper_length > setup->hole_length - setup->hole_bullet_radius ||
      (setup->hole_radius +
       setup->outer_taper_width +
       setup->inner_taper_width) > setup->xtal_radius) {
    printf("\nERROR: Inconsistent detector dimensions:\n"
           "   crystal length and radius: %5.2f %5.2f\n"
           "      hole length and radius: %5.2f %5.2f\n"
           "outer taper length and width: %5.2f %5.2f\n"
           "inner taper length and width: %5.2f %5.2f\n\n",
           setup->xtal_length, setup->xtal_radius,
           setup->hole_length, setup->hole_radius,
           setup->outer_taper_length, setup->outer_taper_width,
           setup->inner_taper_length, setup->inner_taper_width);
    return 1;
  }

  return 0;
}

/* ************************* mjd_fieldgen.h *************************** */

#define MAX_ITS 50000     // default max number of iterations for relaxation

static int report_config(FILE *fp_out, char *config_file_name);
static int grid_init(MJD_Siggen_Setup *setup);
static int ev_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup);
static int wp_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup);
static int write_ev(MJD_Siggen_Setup *setup);
static int write_wp(MJD_Siggen_Setup *setup);
static int do_relax(MJD_Siggen_Setup *setup, int ev_calc);
static int ev_relax_undep(MJD_Siggen_Setup *setup);
static int wp_relax_undep(MJD_Siggen_Setup *setup);
static int interpolate(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup);

/* -------------------------------------- main ------------------- */
int main(char* config_file, int w, int p, int d, float b, char* impurity_profile)
{
//   char **argv;
  // int argc = 6;
  // char* argv[7] = {"asd", "mehdi.config", "-w", "1", "-p", "1", NULL }; /* note that this is possibly  dangerous */
  MJD_Siggen_Setup setup, setup1, setup2;

  float BV;      // bias voltage
  int   WV = 0;  // 0: do not write the V and E values to ppc_ev.dat
                 // 1: write the V and E values to ppc_ev.dat
                 // 2: write the V and E values for both +r, -r (for gnuplot, NOT for siggen)
  int   WD = 0;  // 0: do not write out depletion surface
                 // 1: write out depletion surface to depl_<HV>.dat

  int   i, j;
  FILE  *fp;
  
  printf("Usage: %s <config_file_name> [options]\n"
          "   Possible options:\n"
    "      -b bias_volts\n"
    "      -w {0,1,2}  do_not/do write the field file)\n"
    "      -d {0,1}  do_not/do write the depletion surface)\n"
    "      -p {0,1}  do_not/do write the WP file)\n"
          "      -r rho_spectrum_file_name\n", config_file);
  printf("Possible options:\n"
      "      -b bias_volts\n"
      "      -w {0,1,2} (for WV options)\n"
      "      -p {0,1}   (for WP options)\n"
      "      -r rho_spectrum_file_name\n");
  if (read_config(config_file, &setup)){
    return 1;
  }

  
  strncpy(setup.config_file_name, config_file, sizeof(setup.config_file_name));

  if (setup.xtal_grid < 0.001) setup.xtal_grid = 0.5;
  BV = setup.xtal_HV;
  WV = setup.write_field;
  setup.rho_z_spe[0] = 0;

  if (b != 0) {
    BV = setup.xtal_HV = b;   // bias volts
  }
  if (w >= 0)
    WV = w;               // write-out options
  if (d >= 0)
    WD = d;               // write-out options
  if (p >= 0)
    setup.write_WP = p;   // weighting-potential options
  if (0) {
  // if (impurity_profile) {
    if (!(fp = fopen(impurity_profile, "r"))) {   // impurity-profile-spectrum file name
      printf("\nERROR: cannot open impurity profile spectrum file %s\n\n", impurity_profile);
      return 1;
    }
    fread(setup.rho_z_spe, 36, 1, fp);
    for (j=0; j<1024; j++) setup.rho_z_spe[i] = 0;
    fread(setup.rho_z_spe, sizeof(setup.rho_z_spe), 1, fp);
    fclose(fp);
    printf(" z(mm)   rho\n");
    for (j=0; j < 200 && setup.rho_z_spe[j] != 0.0f; j++)
      printf(" %3d  %7.3f\n", j, setup.rho_z_spe[j]);
  }
  
  if (setup.xtal_length/setup.xtal_grid * setup.xtal_radius/setup.xtal_grid > 2500*2500) {
    printf("Error: Crystal size divided by grid size is too large!\n");
    return 1;
  }
  if (WV < 0 || WV > 2) WV = 0;

  /* -------------- give details of detector geometry */
  if (setup.verbosity >= CHATTY) {
    printf("\n\n"
           "      Crystal: Radius x Length: %.1f x %.1f mm\n",
	   setup.xtal_radius, setup.xtal_length);
    if (setup.hole_length > 0) {
      if (setup.inner_taper_length > 0)
        printf("    Core hole: Radius x length: %.1f x %.1f mm,"
               " taper %.1f x %.1f mm (%2.f degrees)\n",
               setup.hole_radius, setup.hole_length,
               setup.inner_taper_width, setup.inner_taper_length, setup.taper_angle);
      else
        printf("    Core hole: Radius x length: %.1f x %.1f mm\n",
               setup.hole_radius, setup.hole_length);
    }
    printf("Point contact: Radius x length: %.1f x %.1f mm\n",
           setup.pc_radius, setup.pc_length);
    if (setup.ditch_depth > 0) {
      printf("  Wrap-around: Radius x ditch x gap:  %.1f x %.1f x %.1f mm\n",
             setup.wrap_around_radius, setup.ditch_depth, setup.ditch_thickness);
    }
    printf("         Bias: %.0f V\n", BV);
  }
    
  if ((BV < 0 && setup.impurity_z0 < 0) || (BV > 0 && setup.impurity_z0 > 0)) {
    printf("ERROR: Expect bias and impurity to be opposite sign!\n");
    return 1;
  } 
  if (setup.impurity_z0 > 0) {
    // swap polarity for n-type material; this lets me assume all voltages are positive
    BV = -BV;
    setup.xtal_HV *= -1.0;
    setup.impurity_z0         *= -1.0;
    setup.impurity_gradient   *= -1.0;
    setup.impurity_quadratic  *= -1.0;
    setup.impurity_surface    *= -1.0;
    setup.impurity_radial_add *= -1.0;
  }
  /* use an adaptive grid; start out coarse and then refine the grid */
  memcpy(&setup1, &setup, sizeof(setup));
  memcpy(&setup2, &setup, sizeof(setup));
  setup1.xtal_grid *= 9.0;
  setup2.xtal_grid *= 3.0;
  if (grid_init(&setup1) != 0 ||
      grid_init(&setup2) != 0 ||
      grid_init(&setup)  != 0) {
    printf("failed to init field calculations\n");
    return 1;
  }

  /* -------------- calculate electric potential/field */
  if (setup.write_field) {
    setup1.write_field = 0; // no need to save intermediate calculations
    setup2.write_field = 0;
    if (setup.xtal_grid > 0.4) {
      ev_calc(&setup2, NULL);
    } else {
      ev_calc(&setup1, NULL);
      ev_calc(&setup2, &setup1);
    }
    ev_calc(&setup, &setup2);
  }

  /* -------------- calculate weighting potential */
  if (setup.write_WP) {
    setup1.write_WP = 0; // no need to save intermediate calculations
    setup2.write_WP = 0;
    if (setup.xtal_grid > 0.4) {
      wp_calc(&setup2, NULL);
    } else {
      wp_calc(&setup1, NULL);
      wp_calc(&setup2, &setup1);
    }
    wp_calc(&setup, &setup2);
  }

  /* -------------- calculate capacitance
     1/2 * epsilon * integral(E^2) = 1/2 * C * V^2
     so    C = epsilon * integral(E^2) / V^2    V = 1 volt
  */
  double esum, esum2, pi=3.14159, Epsilon=(8.85*16.0/1000.0);  // permittivity of Ge in pF/mm
  float  E_r, E_z;
  float  grid = setup.xtal_grid;
  int    r, z, test;
  int    L  = lrint(setup.xtal_length/grid)+1;
  int    R  = lrint(setup.xtal_radius/grid)+1;
  int    LC = lrint(setup.pc_length/grid)+1;
  int    RC = lrint(setup.pc_radius/grid)+1;

  if (setup.write_WP) {
    esum = esum2 = test = 0;
    for (z=1; z<L; z++) {
      for (r=2; r<R; r++) {
        E_r = setup.eps_dr[z][r]/16.0 * (setup.v[1][z][r] - setup.v[1][z][r+1])/(0.1*grid);
        E_z = setup.eps_dz[z][r]/16.0 * (setup.v[1][z][r] - setup.v[1][z+1][r])/(0.1*grid);
        esum += (E_r*E_r + E_z*E_z) * (double) (r-1);
        if ((r == RC   && z <= LC)   || (r <= RC   && z == LC)   ||
            (r == RC+1 && z <= LC+1) || (r <= RC+1 && z == LC+1)) { // average over two different surfaces
          if (setup.point_type[z+1][r+1] == PC) test = 1;
          esum2 += 0.5 * sqrt(E_r*E_r + E_z*E_z) * (double) (r-1);  // 0.5 since averaging over 2 surfaces
        }
      }
    }
    esum  *= 2.0 * pi * 0.01 * Epsilon * pow(grid, 3.0);
    // Epsilon is in pF/mm
    // 0.01 converts (V/cm)^2 to (V/mm)^2, pow() converts to grid^3 to mm3
    esum2 *= 2.0 * pi * 0.1 * Epsilon * pow(grid, 2.0);
    // 0.1 converts (V/cm) to (V/mm),  grid^2 to  mm2
    printf("  >>  Calculated capacitance at %.0f V: %.3lf pF\n", BV, esum);
    if (!test)
      printf("  >>  Alternative calculation of capacitance: %.3lf pF\n", esum2);
  }

  /* -------------- estimate depletion voltage */
  double min = BV, min2 = BV, dV, dW, testv;
  int    vminr=0, vminz=0;
  int    dz[4] = {1, -1, 0, 0}, dr[4] = {0, 0, 1, -1};
  if (setup.write_WP) {
    if (setup.fully_depleted) {
      // find minimum potential
      for (z=1; z<LC+2; z++) {
        for (r=1; r<RC+2; r++) {
          if (setup.vsave[z][r] > 0 &&
              min > setup.vsave[z][r] / (1.0 - setup.v[1][z][r])) {
            min = setup.vsave[z][r] / (1.0 - setup.v[1][z][r]);
          }
        }
      }
      /* check for bubble depletion / pinch-off by seeing how much the bias
         must be reduced for any pixel to be in a local potential minimum  */
      for (z=LC+2; z<L-2; z++) {
        for (r=1; r<R-2; r++) {
          if (setup.point_type[z][r] == INSIDE && setup.v[1][z][r] > 0.0001) {
            testv = -1;
            for (i=0; i<4; i++) {
              if (r==1 && i==2) break;  // do not check dr for r=1 (=0.0)
              dV = setup.vsave[z+dz[i]][r+dr[i]]  - setup.vsave[z][r];  // potential
              dW = setup.v[1][z+dz[i]][r+dr[i]]   - setup.v[1][z][r];   // WP
              if (dW*grid > 0.00001 && dV < 0 && testv < -dV/dW) testv = -dV/dW;
            }
            if (testv >= 0 && min2 > testv) {
              min2 = testv;
              vminr = r; vminz = z;
            }
          }
        }
      }
      if (min2 < min) {
        printf("Estimated pinch-off voltage = %.0f V\n", BV - min);
        printf(" min2 = %.1f at (r,z) = (%.1f, %.1f), so\n",
               min2, (vminr-1)*grid, (vminz-1)*grid);
        printf("   Full depletion (max pinch-off voltage) = %.0f\n", BV - min2);
      } else {
        printf("Estimated depletion voltage = %.0f V\n", BV - min);
      }
    }

    printf("Minimum bulk field = %.2f V/cm at (r,z) = (%.1f, %.1f) mm\n\n",
           setup.Emin, setup.rmin, setup.zmin);
  }
 
  return 0;
} /* main */

/* -------------------------------------- ev_calc ------------------- */
int ev_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup) {
  int    i, j;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+3;
  int    R  = lrint(setup->xtal_radius/grid)+3;

  if (!old_setup) {
    printf("\n\n ---- starting EV calculation --- \n");
    for (i = 1; i < L; i++) {
      for (j = 1; j < R; j++) {
        setup->v[0][i][j] = setup->v[1][i][j] = setup->xtal_HV/2.0;
      }
    }
  }
  if (old_setup) interpolate(setup, old_setup);
  setup->fully_depleted = 1;
  setup->bubble_volts = 0;

  /* set boundary voltages */
  for (i = 1; i < L; i++) {
    for (j = 1; j < R; j++) {
      if (setup->point_type[i][j] == HVC)
        setup->v[0][i][j] = setup->v[1][i][j] = setup->xtal_HV;
      if (setup->point_type[i][j] == PC)
        setup->v[0][i][j] = setup->v[1][i][j] = 0.0;
    }
  }

  if (!old_setup || !old_setup->fully_depleted) ev_relax_undep(setup);
  else do_relax(setup, 1);
  if (setup->write_field) write_ev(setup);

  if (setup->fully_depleted) {
    printf("Detector is fully depleted.\n");
    /* save potential close to point contact, to use later when calculating depletion voltage */
    for (i = 1; i < L; i++) {
      for (j = 1; j < R; j++) {
        setup->vsave[i][j] = fabs(setup->v[1][i][j]);
      }
    }
  } else {
    printf("Detector is not fully depleted.\n");
    if (setup->bubble_volts > 0)
      printf("Pinch-off bubble at %.1f V potential\n", setup->bubble_volts);
    if (!old_setup) {
      // write a little file that shows any undepleted voxels in the crystal
      FILE *file = fopen("undepleted.txt", "w");
      for (j = R-1; j > 0; j--) {
	setup->undepleted[j][L-1] = '\0';
	fprintf(file, "%s\n", setup->undepleted[j]+1);
      }
      fclose(file);
    }
  }

  return 0;
} /* ev_calc */

/* -------------------------------------- wp_calc ------------------- */
int wp_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup) {
  int    i, j, pinched_off = 0;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+3;
  int    R  = lrint(setup->xtal_radius/grid)+3;

  if (!old_setup) {
    printf("\n\n ---- starting WP calculation --- \n");
    for (i = 1; i < L; i++) {
      for (j = 1; j < R; j++) {
        setup->v[0][i][j] = setup->v[1][i][j] = 0.01;
      }
    }
  }
  if (old_setup) interpolate(setup, old_setup);

  /* set boundary voltages */
  for (i = 1; i < L; i++) {
    for (j = 1; j < R; j++) {
      if (setup->point_type[i][j] == HVC)
        setup->v[0][i][j] = setup->v[1][i][j] = 0.0;
      else if (setup->point_type[i][j] == PC)
        setup->v[0][i][j] = setup->v[1][i][j] = 1.0;
      else if (setup->undepleted[j][i] == '*') {
        setup->point_type[i][j] = PC;
        setup->v[0][i][j] = setup->v[1][i][j] = 1.0;
      } else if (setup->undepleted[j][i] == 'B') {
        setup->point_type[i][j] = PINCHOFF;
        pinched_off = 1;
      }
    }
  }

  if (pinched_off) wp_relax_undep(setup);
  else do_relax(setup, 0);

  if (setup->write_WP) write_wp(setup);

  return 0;
} /* wp_calc */

/* -------------------------------------- report_config ------------------- */
int report_config(FILE *fp_out, char *config_file_name) {
  char  *c, line[256];
  FILE  *file;

  fprintf(fp_out, "# Config file: %s\n", config_file_name);
  if (!(file = fopen(config_file_name, "r"))) return 1;

  while (fgets(line, sizeof(line), file)) {
    if (strlen(line) < 3 || *line == ' ' || *line == '\t' || *line == '#') continue;
    if ((c = strchr(line, '#')) || (c = strchr(line, '\n'))) *c = '\0';
    fprintf(fp_out, "# %s\n", line);
  }
  fclose(file);
  return 0;
} /* report_config */

/* -------------------------------------- write_ev ------------------- */
int write_ev(MJD_Siggen_Setup *setup) {
  int    i, j;
  int neww = 1;
  float  grid = setup->xtal_grid;
  //int    L  = setup->xtal_length/grid + 2.99;
  //int    R  = setup->xtal_radius/grid + 2.99;
  int    L  = lrint(setup->xtal_length/grid)+2;
  int    R  = lrint(setup->xtal_radius/grid)+2;
  float  r, z;
  float  E_r, E_z, E;
  FILE   *file;
  double ***v = setup->v;

  setup->Emin = 99999.9;
  setup->rmin = setup->zmin = 999.9;


  if (setup->impurity_z0 > 0) {
    // swap voltages back to negative for n-type material
    for (i=1; i<L; i++) {
      for (j=1; j<R; j++) {
        setup->v[neww][i][j] = -setup->v[neww][i][j];
      }
    }
  }

  /* write potential and field to output file */
  if (!(file = fopen(setup->field_name, "w"))) {
    printf("ERROR: Cannot open file %s for electric field...\n", setup->field_name);
    return 1;
  } else {
    printf("Writing electric field data to file %s\n", setup->field_name);
  }
  /* copy configuration parameters to output file */
  report_config(file, setup->config_file_name);
  fprintf(file, "#\n# HV bias in fieldgen: %.1f V\n", setup->xtal_HV);
  if (setup->fully_depleted) {
    fprintf(file, "# Detector is fully depleted.\n");
  } else {
    fprintf(file, "# Detector is not fully depleted.\n");
    if (setup->bubble_volts > 0.0f)
      fprintf(file, "# Pinch-off bubble at %.0f V potential\n", setup->bubble_volts);
  }
  fprintf(file, "#\n## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)\n");
  
  for (j = 1; j < R; j++) {
    r = (j-1) * grid;
    for (i = 1; i < L; i++) {
      z = (i-1) * grid;
      // calc E in r-direction
      if (j == 1) {  // r = 0; symmetry implies E_r = 0
        E_r = 0;
      } else if (setup->point_type[i][j] == CONTACT_EDGE) {
        E_r = ((v[neww][i][j] - v[neww][i][j+1])*setup->dr[1][i][j] +
               (v[neww][i][j-1] - v[neww][i][j])*setup->dr[0][i][j]) / (0.2*grid);
      } else if (setup->point_type[i][j] < INSIDE &&
                 setup->point_type[i][j-1] == CONTACT_EDGE) {
        E_r =  (v[neww][i][j-1] - v[neww][i][j]) * setup->dr[1][i][j-1] / ( 0.1*grid) ;
      } else if (setup->point_type[i][j] < INSIDE &&
                 setup->point_type[i][j+1] == CONTACT_EDGE) {
        E_r =  (v[neww][i][j] - v[neww][i][j+1]) * setup->dr[0][i][j+1] / ( 0.1*grid) ;
      } else if (j == R-1) {
        E_r = (v[neww][i][j-1] - v[neww][i][j])/(0.1*grid);
      } else {
        E_r = (v[neww][i][j-1] - v[neww][i][j+1])/(0.2*grid);
      }
      // calc E in z-direction
      if (setup->point_type[i][j] == CONTACT_EDGE) {
        E_z = ((v[neww][i][j] - v[neww][i+1][j])*setup->dz[1][i][j] +
               (v[neww][i-1][j] - v[neww][i][j])*setup->dz[0][i][j]) / (0.2*grid);
      } else if (setup->point_type[i][j] < INSIDE &&
                 setup->point_type[i-1][j] == CONTACT_EDGE) {
        E_z =  (v[neww][i-1][j] - v[neww][i][j]) * setup->dz[1][i-1][j] / ( 0.1*grid) ;
      } else if (setup->point_type[i][j] < INSIDE &&
                 setup->point_type[i+1][j] == CONTACT_EDGE) {
        E_z =  (v[neww][i][j] - v[neww][i+1][j]) * setup->dz[0][i+1][j] / ( 0.1*grid) ;
      } else if (i == 1) {
        E_z = (v[neww][i][j] - v[neww][i+1][j])/(0.1*grid);
      } else if (i == L-1) {
        E_z = (v[neww][i-1][j] - v[neww][i][j])/(0.1*grid);
      } else {
        E_z = (v[neww][i-1][j] - v[neww][i+1][j])/(0.2*grid);
      }
      E = sqrt(E_r*E_r + E_z*E_z);
      fprintf(file, "%7.2f %7.2f %7.1f %7.1f %7.1f %7.1f\n",
              r, z, v[neww][i][j], E, E_r, E_z);

      /* check for minimum field inside bulk of detector */
      int k = 3.0/grid;
      if (E > 0.1 && E < setup->Emin &&
          i > k+1 && j < R-k-1 && i < L-k-1 &&
          setup->point_type[i][j] == INSIDE &&
          setup->point_type[i + k][j] == INSIDE &&  // point is at least 3 mm from a boundary
          setup->point_type[i - k][j] == INSIDE &&
          setup->point_type[i][j + k] == INSIDE &&
          (j < k+1 || setup->point_type[i][j - k] == INSIDE)) {
        setup->Emin = E;
        setup->rmin = r;
        setup->zmin = z;
      }
    }
    fprintf(file, "\n");
  }
  fclose(file);
  if (!setup->write_WP)
    printf("\n Minimum bulk field = %.2f V/cm at (r,z) = (%.1f, %.1f) mm\n\n",
           setup->Emin, setup->rmin, setup->zmin);

  if (0) { /* write point_type to output file */
    file = fopen("fields/point_type.dat", "w");
    for (j = 1; j < R; j++) {
      for (i = 1; i < L; i++)
        fprintf(file, "%7.2f %7.2f %2d\n",
                (j-1)*grid, (i-1)*grid, setup->point_type[i][j]);
      fprintf(file, "\n");
    }
    fclose(file);
  }

  return 0;
 } /* write_ev */

/* -------------------------------------- write_wp ------------------- */
 int write_wp(MJD_Siggen_Setup *setup) {
  int    i, j, neww=1;;
  float  grid = setup->xtal_grid;
  int    L  = setup->xtal_length/grid + 2.99;
  int    R  = setup->xtal_radius/grid + 2.99;
  //int    L  = lrint(setup->xtal_length/grid)+2;
  //int    R  = lrint(setup->xtal_radius/grid)+2;
  float  r, z;
  FILE *file;

  if (!(file = fopen(setup->wp_name, "w"))) {
    printf("ERROR: Cannot open file %s for weighting potential...\n", setup->wp_name);
    return 1;
  } else {
    printf("Writing weighting potential to file %s\n\n", setup->wp_name);
  }

  /* copy configuration parameters to output file */
  report_config(file, setup->config_file_name);
  fprintf(file, "#\n# HV bias in fieldgen: %.1f V\n", setup->xtal_HV);
  if (setup->fully_depleted) {
    fprintf(file, "# Detector is fully depleted.\n");
  } else {
    fprintf(file, "# Detector is not fully depleted.\n");
    if (setup->bubble_volts > 0.0f)
      fprintf(file, "# Pinch-off bubble at %.0f V potential\n", setup->bubble_volts);
  }
  fprintf(file, "#\n## r (mm), z (mm), WP\n");
  for (j = 1; j < R; j++) {
    r = (j-1) * grid;
    for (i = 1; i < L; i++) {
      z = (i-1) * grid;
      fprintf(file, "%7.2f %7.2f %12.6e\n", r, z, setup->v[neww][i][j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  return 0;
 } /* write_wp */

/* -------------------------------------- dist_from_contact ------------------- */
float dist_from_contact(cyl_pt pt, cyl_pt delta, MJD_Siggen_Setup *setup) {
  float  factor = 1, d = 0.5;
  cyl_pt test;
  int    n;

  for (n=0; n<7; n++) {  // 7 steps => 1/128 precision
    test.r = pt.r + factor * delta.r;
    test.z = pt.z + factor * delta.z;
    if (outside_detector_cyl(test, setup)) {
      factor -= d;
    } else {
      if (n == 0) return -1.0;
      factor += d;
    } 
    d /= 2.0;
  }
  return factor;
} /* dist_from_contact */

/* -------------------------------------- grid_init ------------------- */
#define SQ(x) ((x)*(x))
int grid_init(MJD_Siggen_Setup *setup) {
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+3;
  int    R  = lrint(setup->xtal_radius/grid)+3;
  int    i, j;
  float  r, z;


  /* first malloc arrays in setup */
  if ((setup->impurity   = (double**)malloc(L * sizeof(*setup->impurity)))   == NULL ||
      (setup->eps        = (double**)malloc(L * sizeof(*setup->eps)))        == NULL ||
      (setup->eps_dr     = (double**)malloc(L * sizeof(*setup->eps_dr)))     == NULL ||
      (setup->eps_dz     = (double**)malloc(L * sizeof(*setup->eps_dz)))     == NULL ||
      (setup->v[0]       = (double**)malloc(L * sizeof(*setup->v[0])))       == NULL ||
      (setup->v[1]       = (double**)malloc(L * sizeof(*setup->v[1])))       == NULL ||
      (setup->dr[0]      = (double**)malloc(L * sizeof(*setup->dr[0])))      == NULL ||
      (setup->dr[1]      = (double**)malloc(L * sizeof(*setup->dr[1])))      == NULL ||
      (setup->dz[0]      = (double**)malloc(L * sizeof(*setup->dz[0])))      == NULL ||
      (setup->dz[1]      = (double**)malloc(L * sizeof(*setup->dz[1])))      == NULL ||
      (setup->undepleted = (char**)malloc(R * sizeof(*setup->undepleted))) == NULL ||
      (setup->s1         = (double*)malloc(R * sizeof(*setup->s1)))         == NULL ||
      (setup->s2         = (double*)malloc(R * sizeof(*setup->s2)))         == NULL ||
      (setup->vsave      = (double**)malloc(L * sizeof(*setup->vsave)))      == NULL ||
      (setup->point_type = (char**)malloc(L * sizeof(*setup->point_type))) == NULL) {
    printf("malloc failed\n");
    return -1;
  }
  /* start from i=1 so that i=0 can be used for reflection symmetry around r=0 or z=0 */
  for (i = 1; i < L; i++) {
    if ((setup->impurity[i]   = (double*)malloc(R * sizeof(**setup->impurity)))   == NULL ||
        (setup->v[0][i]       = (double*)malloc(R * sizeof(**setup->v[0])))       == NULL ||
	      (setup->v[1][i]       = (double*)malloc(R * sizeof(**setup->v[1])))       == NULL ||
        (setup->dr[0][i]      = (double*)malloc(R * sizeof(**setup->dr[0])))      == NULL ||
	      (setup->dr[1][i]      = (double*)malloc(R * sizeof(**setup->dr[1])))      == NULL ||
        (setup->dz[0][i]      = (double*)malloc(R * sizeof(**setup->dz[0])))      == NULL ||
        (setup->dz[1][i]      = (double*)malloc(R * sizeof(**setup->dz[1])))      == NULL ||
        (setup->vsave[i]      = (double*)malloc(R * sizeof(**setup->vsave)))      == NULL ||
        (setup->point_type[i] = (char*)malloc(R * sizeof(**setup->point_type))) == NULL) {
      printf("malloc failed\n");
      return -1;
    }
  }
  for (i = 1; i < L; i++) {
    if ((setup->eps[i]        = (double*)malloc(R * sizeof(**setup->eps)))    == NULL ||
        (setup->eps_dr[i]     = (double*)malloc(R * sizeof(**setup->eps_dr))) == NULL ||
        (setup->eps_dz[i]     = (double*)malloc(R * sizeof(**setup->eps_dz))) == NULL) {
      printf("malloc failed\n");
      return -1;
    }
    for (j = 0; j < R; j++)
      setup->eps[i][j] = setup->eps_dz[i][j] = setup->eps_dr[i][j] = 16.0;
  }
  for (j = 0; j < R; j++) {
    if ((setup->undepleted[j] = (char*)malloc(L * sizeof(**setup->undepleted))) == NULL) {
      printf("malloc failed\n");
      return -1;
    }
    memset(setup->undepleted[j], ' ', L);
  }

  /* set up reflection symmetry around r=0 or z=0 */
  setup->impurity[0]   = (double*)malloc(R * sizeof(**setup->impurity));
  setup->v[0][0]   = setup->v[0][2];
  setup->v[1][0]   = setup->v[1][2];
  setup->eps[0]    = setup->eps[1];
  setup->eps_dr[0] = setup->eps_dr[1];
  setup->eps_dz[0] = setup->eps_dz[1];
  setup->point_type[0] = setup->point_type[1];

  /* ------------------------------------------------------------ */
  /* weighting values for the relaxation alg. as a function of r
     in the following we divide areas and volumes by pi
     r_bin   rmax  A_top A_outside A_inside  volume  total_surf  out/top  tot/vol
     0     1/2    1/4      1         0       1/4      1.5         4        6  << special case
     1     3/2      2      3         1        2        8        3/2        4
     2     5/2      4      5         3        4       16        5/4        4
     3     7/2      6      7         5        6       24        7/6        4
     r   r+0.5     2r    2r+1      2r-1      2r       8r     (2r+1)/2r     4
     = 1+0.5/r
  */
  setup->s1[1] = 4.0;
  setup->s2[1] = 0.0;
  for (i=2; i<R; i++) {
    setup->s1[i] = 1.0 + 0.5 / (double) (i-1);   //  for r+1
    setup->s2[i] = 1.0 - 0.5 / (double) (i-1);   //  for r-1
  }
  setup->s2[1] = setup->s1[1]; // special case for reflection symm at r=0

  for (i = 1; i < L; i++) {
    for (j = 0; j < R; j++) {
      setup->dr[0][i][j] = setup->s2[j];   //  for r-1
      setup->dr[1][i][j] = setup->s1[j];   //  for r+1
      setup->dz[0][i][j] = 1;       //  for z-1
      setup->dz[1][i][j] = 1;       //  for z+1
    }
  }

  /* set up pixel point types for boundary conditions etc */
  float  d, g = 0.05 * grid, lith  = setup->Li_thickness;
  cyl_pt pt, pt1, pt2, pt3, pt4;

  setup->rmax = setup->xtal_radius - lith;
  setup->zmax = setup->xtal_length - lith;
  setup->hole_radius += lith;
  setup->hole_bullet_radius += lith;
  setup->bottom_taper_length += 0.71*lith; // add top tapers?

  for (i = 1; i < L; i++) {
    pt.z = pt3.z = pt4.z = (i-1) * grid;
    pt1.z = pt.z + g;
    pt2.z = pt.z - g;
    for (j = 1; j < R; j++) {
      pt.r = pt1.r = pt2.r = (j-1) * grid;
      pt3.r = pt.r + g;
      pt4.r = pt.r - g;
      setup->point_type[i][j] = INSIDE;

      /* see if pixel is ouside (or very nearly outside) the detector bulk */
      if (i == 1 || (pt.r >= setup->wrap_around_radius && pt.z - g < lith) ||
          outside_detector_cyl(pt1, setup) || outside_detector_cyl(pt2, setup) ||
          outside_detector_cyl(pt3, setup) || outside_detector_cyl(pt4, setup)) {

        /* check for inside ditch
           boundary condition at Ge-vacuum interface:
           epsilon0 * E_vac = espilon_Ge * E_Ge  */
        if (setup->ditch_depth > 0 && pt.z < setup->ditch_depth + grid &&
            pt.r <= setup->wrap_around_radius &&
            pt.r >= setup->wrap_around_radius - setup->ditch_thickness) {
          setup->point_type[i][j] = DITCH;
          setup->eps[i][j] = setup->eps_dz[i][j] = setup->eps_dr[i][j] = 1.0;

          /* check for inside (point) contact */
        } else if (pt.z < setup->pc_length + g && pt.r < setup->pc_radius+g) {
          setup->point_type[i][j] = PC;

        /* check for passivated area  */
        } else if (i == 1 && pt.r < setup->wrap_around_radius &&      // BEGE/ICPC
                   pt.r < setup->rmax - setup->bottom_taper_length) { // PPC
          setup->point_type[i][j] = PASSIVE;

        /* only remaining surface is HV contact */
        } else  {
          setup->point_type[i][j] = HVC;
        }
      }
    }
  }

  /* find the pixels next to the contact surfaces */
  cyl_pt dp1, dp2, dp3, dp4;
  dp1.z = dp3.r = grid;
  dp2.z = dp4.r = -grid;
  dp1.r = dp2.r = dp3.z = dp4.z = 0;
  for (i = 2; i < L-1; i++) {
    pt.z = (i-1) * grid;
    for (j = 1; j < R; j++) {
      pt.r = (j-1) * grid;
      if (setup->point_type[i][j] == INSIDE &&
          (setup->point_type[i+1][j] < INSIDE || setup->point_type[i-1][j] < INSIDE ||
           setup->point_type[i][j+1] < INSIDE || (j > 1 && setup->point_type[i][j-1] < INSIDE))) {
        setup->point_type[i][j] = CONTACT_EDGE;
        /* find distance to contact surface */
        if (setup->point_type[i+1][j] < INSIDE && (d = dist_from_contact(pt, dp1, setup)) > 0)
          setup->dz[1][i][j] = 1.0/d;
        if (setup->point_type[i-1][j] < INSIDE && (d = dist_from_contact(pt, dp2, setup)) > 0)
          setup->dz[0][i][j] = 1.0/d;
        else if (setup->point_type[i-1][j] < INSIDE &&
                 pt.r >= setup->wrap_around_radius && pt.z - grid < lith)
          setup->dz[0][i][j] = grid/(pt.z - lith);
        if (setup->point_type[i][j+1] < INSIDE && (d = dist_from_contact(pt, dp3, setup)) > 0)
          setup->dr[1][i][j] = setup->s1[j] * 1.0/d;
        if (j > 1 && setup->point_type[i][j-1] < INSIDE &&
            (d = dist_from_contact(pt, dp4, setup)) > 0)
          setup->dr[0][i][j] = setup->s2[j] * 1.0/d;
      }
    }
  }
  setup->hole_radius -= lith;
  setup->hole_bullet_radius -= lith;
  setup->bottom_taper_length -= 0.71*lith;

  /* for pixels adjacent to the ditch, set point_type to DITCH_EDGE
     and for z=0, set flag for passivated surface */
  for (i = 1; i < L; i++) {
    for (j = 1; j < R; j++) {
      setup->eps_dr[i][j-1] = (setup->eps[i][j-1] + setup->eps[i][j]) / 2.0f;
      setup->eps_dz[i-1][j] = (setup->eps[i-1][j] + setup->eps[i][j]) / 2.0f;
      if (setup->point_type[i][j] == INSIDE &&
          (setup->point_type[i-1][j] == DITCH ||
           setup->point_type[i][j-1] == DITCH ||
           setup->point_type[i][j+1] == DITCH)) setup->point_type[i][j] = DITCH_EDGE;
    }
    setup->eps_dr[i][0] = setup->eps_dr[i][1];
  }

  /* set up impurity array */
  double *imp_z, imp_ra = 0, imp_rm = 1;
  double e_over_E = 11.310; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0
  /* malloc local array */
  if ((imp_z  = (double*)malloc(L * sizeof(*imp_z))) == NULL) {
    printf("malloc failed\n");
    return -1;
  }

  if (setup->rho_z_spe[0] == 0) {
    for (i = 1; i < L; i++) {
      z = (i-1) * grid;
      imp_z[i] = e_over_E * grid*grid / 4.0 *
                 (setup->impurity_z0 +
                  setup->impurity_gradient * z * 0.1 +
                  setup->impurity_quadratic *
                  (1.0 - SQ(z - setup->xtal_length/2.0) / SQ(setup->xtal_length/2.0)));
    }
  } else {
    for (i = 1; i < L; i++)
      imp_z[i] = e_over_E * grid*grid / 4.0 * setup->rho_z_spe[(int) ((i-1) * grid)];
  }
  for (j = 1; j < R; j++) {
    r = (j-1) * grid;
    if (setup->impurity_rpower > 0.1) {
      imp_ra = setup->impurity_radial_add * e_over_E * grid*grid / 4.0 *
        pow((double) r / setup->xtal_radius, setup->impurity_rpower);
      imp_rm = 1.0 + (setup->impurity_radial_mult - 1.0f) *
        pow((double) r / setup->xtal_radius, setup->impurity_rpower);
    }
    for (i = 1; i < L; i++)  setup->impurity[i][j] = imp_z[i] * imp_rm + imp_ra;
    if (setup->point_type[1][j] == PASSIVE) {
      setup->impurity[1][j] += setup->impurity_surface * e_over_E * grid/4.0;
    }
    /* reduce charge volume for CONTACT_EDGE pixels */
    for (i = 1; i < L; i++) {
      if (setup->point_type[i][j] == CONTACT_EDGE) {
        setup->impurity[i][j] /=
          SQ(setup->dz[1][i][j]*setup->dz[0][i][j] * setup->dr[1][i][j]*setup->dr[0][i][j]);
      }
    }
  }
           
  /* free local and no-longer-needed arrays */
  free(imp_z);
  for (i = 1; i < L; i++) free(setup->eps[i]);
  free(setup->eps);

  return 0;   
} /* grid_init */
#undef SQ

/* -------------------------------------- do_relax ------------------- */
int do_relax(MJD_Siggen_Setup *setup, int ev_calc) {
  int    old = 1, neww = 0, iter, r, z;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+2;
  int    R  = lrint(setup->xtal_radius/grid)+2;
  double eps_sum, v_sum, dif, sum_dif, max_dif;
  double ***v = setup->v, **eps_dr = setup->eps_dr, **eps_dz = setup->eps_dz;
  double ***dr = setup->dr, ***dz = setup->dz;
  double *s1 = setup->s1, *s2 = setup->s2;
  double e_over_E = 11.310; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0


  if (ev_calc) {
    // for field calculation, save impurity value along passivated surface
    for (r = 1; r < R; r++)
      setup->impurity[0][r] = setup->impurity[1][r];
  } else {
    // for WP calculation, clear all impurity values
    for (z = 0; z < L; z++) {
      for (r = 1; r < R; r++) {
        setup->impurity[z][r] = 0;
      }
    }
  }

  for (iter = 0; iter < setup->max_iterations; iter++) {

    /* the following definition of the factor for over-relaxation improves convergence
           time by a factor ~ 70-120 for a 2kg ICPC detector, grid = 0.1 mm
         OR_fact increases with increasing volxel count (L*R)
               and with increasing iteration number
         0.997 is maximum asymptote for very large pixel count and iteration number */
    double OR_fact;
    if (ev_calc)  OR_fact = (1.991 - 1500.0/(L*R));
    else          OR_fact = (1.992 - 1500.0/(L*R));
    if (OR_fact < 1.4) OR_fact = 1.4;
    // if (iter == 0) printf("OR_fact = %f\n", OR_fact);
    if (iter < 1) OR_fact = 1.0;

    old = neww;
    neww = 1 - neww;
    sum_dif = 0;
    max_dif = 0;

    if (setup->vacuum_gap > 0) {   // modify impurity value along passivated surface
      for (r = 1; r < R; r++)      //   due to surface charge induced by capacitance
        setup->impurity[1][r] = setup->impurity[0][r] -
          v[old][1][r] * 5.52e-4 * e_over_E * grid / setup->vacuum_gap;
    }

    /* start from z=1 and r=1 so that (z,r)=0 can be
       used for reflection symmetry around r=0 or z=0 */
    for (z = 1; z < L; z++) {
      /* manage r=0 reflection symmetry */
      setup->v[old][z][0] = setup->v[neww][z][0] = setup->v[old][z][2];

      for (r = 1; r < R; r++) {
        if (setup->point_type[z][r] < INSIDE) continue;   // HV or point contact

        if (setup->point_type[z][r] < DITCH) {       // normal bulk or passivated surface, no complications
          v_sum = (v[old][z+1][r] + v[old][z][r+1]*s1[r] +
                   v[neww][z-1][r] + v[neww][z][r-1]*s2[r]);
          if (r > 1) eps_sum = 4;
          else       eps_sum = 2 + s1[r] + s2[r];
        } else if (setup->point_type[z][r] == CONTACT_EDGE) {  // adjacent to the contact
          v_sum = (v[old][z+1][r]*dz[1][z][r] + v[old][z][r+1]*dr[1][z][r] +
                   v[neww][z-1][r]*dz[0][z][r] + v[neww][z][r-1]*dr[0][z][r]);
          eps_sum = dz[1][z][r] + dr[1][z][r] + dz[0][z][r] + dr[0][z][r];
        } else if (setup->point_type[z][r] >= DITCH) {  // in or adjacent to the ditch
          v_sum = (v[old][z+1][r]*eps_dz[z  ][r] + v[old][z][r+1]*eps_dr[z][r  ]*s1[r] +
                   v[neww][z-1][r]*eps_dz[z-1][r] + v[neww][z][r-1]*eps_dr[z][r-1]*s2[r]);
          eps_sum = (eps_dz[z][r]   + eps_dr[z][r]  *s1[r] +
                     eps_dz[z-1][r] + eps_dr[z][r-1]*s2[r]);
        }

        // calculate the interpolated mean potential and the effect of the space charge

        if ((ev_calc || (setup->vacuum_gap > 0 && z == 1)) &&
            setup->point_type[z][r] < CONTACT_EDGE && r > 1 && z > 1) {   // normal bulk, no complications
          v[neww][z][r] = (1.0-OR_fact)*v[old][z][r] + OR_fact * (v_sum / eps_sum + setup->impurity[z][r]);
        } else if (ev_calc || (setup->vacuum_gap > 0 && z == 1)) {
          v[neww][z][r] = v_sum / eps_sum + setup->impurity[z][r];
        } else if (setup->point_type[z][r] < CONTACT_EDGE && r > 1 && z > 1) {   // normal bulk, no complications
          v[neww][z][r] = (1.0-OR_fact)*v[old][z][r] + OR_fact * v_sum / eps_sum;
        } else {                          // over-relaxation at the edges seems to make things worse
          v[neww][z][r] = v_sum / eps_sum;
        }

        // calculate difference from last iteration, for convergence check
        dif = fabs(v[old][z][r] - v[neww][z][r]);
        sum_dif += dif;
        if (max_dif < dif) max_dif = dif;
      }
    }

    // report results for some iterations
    if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0) {
      if (0 && ev_calc) {
        printf("%5d %d %d %.10f %.10f\n", iter, old, neww, max_dif, sum_dif/(L-2)/(R-2));
      } else {
        printf("%5d %d %d %.10f %.10f ; %.10f %.10f\n",
               iter, old, neww, max_dif, sum_dif/(L-2)/(R-2),
               v[neww][L/2][R/2], v[neww][L/3][R/3]);
      }
    }
    // check for convergence
    if ( ev_calc && max_dif < 0.00000008) break;
    if ( ev_calc && max_dif < 0.0008) break;  // comment out if you want convergence at the numerical error level
    if (!ev_calc && max_dif < 0.0000000001) break;
    if (!ev_calc && max_dif < 0.000001) break;  // comment out if you want convergence at the numerical error level

    /* every 100 iterations, check that detector is really depleted*/
    if (ev_calc && iter > 190 && iter%100 == 0) {
      for (z = 1; z < L; z++) {
        setup->v[old][z][0] = setup->v[old][z][2];
        for (r = 1; r < R; r++) {
          if (setup->point_type[z][r] < INSIDE) continue;   // HV or point contact
          if (v[neww][z][r] < 0 ||
              (v[neww][z][r] < v[neww][z][r] &&
               v[neww][z][r] < v[neww][z][r] &&
               v[neww][z][r] < v[neww][z][r] &&
               v[neww][z][r] < v[neww][z][r])) {
            printf("Detector may not be fully depleted. Switching to ev_relax_undep()\n");
            ev_relax_undep(setup);
            return 0;
          }
        }
      }

    }
  }

  printf(">> %d %.16f\n\n", iter, sum_dif);
  if (setup->vacuum_gap > 0) {   // restore impurity value along passivated surface
    for (r = 1; r < R; r++)
      setup->impurity[1][r] = setup->impurity[0][r];
  }

  return 0;
} /* do_relax */

/* -------------------------------------- ev_relax_undep ------------------- */
/*  This function, unlike do_relax() above, properly handles undepleted detectors.
    Note that this function uses a modified sequential over-relaxtion algorithm,
    while do_relax() above uses a more standard text-book version.
 */
int ev_relax_undep(MJD_Siggen_Setup *setup) {
  int    old = 1, neww = 0, iter, r, z, bvn;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+2;
  int    R  = lrint(setup->xtal_radius/grid)+2;
  double eps_sum, v_sum, save_dif, min;
  double dif, sum_dif, max_dif, bubble_volts;
  double ***v = setup->v, **eps_dr = setup->eps_dr, **eps_dz = setup->eps_dz;
  double ***dr = setup->dr, ***dz = setup->dz;
  double *s1 = setup->s1, *s2 = setup->s2;
  char   **undep = setup->undepleted;
  double e_over_E = 11.310; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0


  // save impurity value along passivated surface
  for (r = 1; r < R; r++)
    setup->impurity[0][r] = setup->impurity[1][r];

  /* initialise the undepleted array for use with bubble depletion */
  for (z = 1; z < L; z++) {
    for (r = 1; r < R; r++) {
      if (setup->point_type[z][r] >= INSIDE) undep[r][z] = 0;
    }
  }

  for (iter = 0; iter < setup->max_iterations; iter++) {

    double OR_fact = ((0.997 - 300.0/(L*R)) * (1.0 - 0.9/(double)(1+iter/6)));
    if (300.0/(L*R) > 0.5) OR_fact = (0.5 * (1.0 - 0.9/(double)(1+iter/6)));
    if (iter < 2) OR_fact = 0.0;

    old = neww;
    neww = 1 - neww;
    sum_dif = 0;
    max_dif = 0;
    bubble_volts = 0;
    bvn = 0;

    if (setup->vacuum_gap > 0) {   // modify impurity value along passivated surface
      for (r = 1; r < R; r++)      //   due to surface charge induced by capacitance
        setup->impurity[1][r] = setup->impurity[0][r] -
          v[old][1][r] * 5.52e-4 * e_over_E * grid / setup->vacuum_gap;
    }

    /* start from z=1 and r=1 so that (z,r)=0 can be
       used for reflection symmetry around r=0 or z=0 */
    for (z = 1; z < L; z++) {
      /* manage r=0 reflection symmetry */
      setup->v[old][z][0] = setup->v[old][z][2];

      for (r = 1; r < R; r++) {
        if (setup->point_type[z][r] < INSIDE) continue;   // HV or point contact
        save_dif = v[old][z][r] - v[neww][z][r];      // step difference from previous iteration

        if (setup->point_type[z][r] < DITCH) {       // normal bulk or passivated surface, no complications
          v_sum = (v[old][z+1][r] + v[old][z][r+1]*s1[r] +
                   v[old][z-1][r] + v[old][z][r-1]*s2[r]);
          if (r > 1) eps_sum = 4;
          else       eps_sum = 2 + s1[r] + s2[r];
        } else if (setup->point_type[z][r] == CONTACT_EDGE) {  // adjacent to the contact
          v_sum = (v[old][z+1][r]*dz[1][z][r] + v[old][z][r+1]*dr[1][z][r] +
                   v[old][z-1][r]*dz[0][z][r] + v[old][z][r-1]*dr[0][z][r]);
          eps_sum = dz[1][z][r] + dr[1][z][r] + dz[0][z][r] + dr[0][z][r];
        } else if (setup->point_type[z][r] >= DITCH) {  // in or adjacent to the ditch
          v_sum = (v[old][z+1][r]*eps_dz[z  ][r] + v[old][z][r+1]*eps_dr[z][r  ]*s1[r] +
                   v[old][z-1][r]*eps_dz[z-1][r] + v[old][z][r-1]*eps_dr[z][r-1]*s2[r]);
          eps_sum = (eps_dz[z][r]   + eps_dr[z][r]  *s1[r] +
                     eps_dz[z-1][r] + eps_dr[z][r-1]*s2[r]);
        }

        // calculate the interpolated mean potential and the effect of the space charge
        min = fminf(fminf(v[old][z+1][r], v[old][z][r+1]),
                    fminf(v[old][z-1][r], v[old][z][r-1]));
        v[neww][z][r] = v_sum / eps_sum + setup->impurity[z][r];

        undep[r][z] /= 2;
        if (v[neww][z][r] <= 0) {
          v[neww][z][r] = 0;
          undep[r][z] = 4;  // do not do over-relaxation for 3 iterations
        } else if (v[neww][z][r] <= min) {
          if (bubble_volts == 0) bubble_volts = min + 0.2*grid*grid; // finer grids require smaller increment here
          v[neww][z][r] = bubble_volts;
          bvn++;
          undep[r][z] = 8;  // do not do over-relaxation for 4 iterations
        }

        // calculate difference from last iteration, for convergence check
        dif = v[old][z][r] - v[neww][z][r];
        if (dif < 0) dif = -dif;
        sum_dif += dif;
        if (max_dif < dif) max_dif = dif;
        // do over-relaxation
        if (!undep[r][z])  v[neww][z][r] += OR_fact*save_dif;
      }
    }

    // report results for some iterations
    if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0) {
      if (0) {
        printf("%5d %d %d %.10f %.10f\n", iter, old, neww, max_dif, sum_dif/(L-2)/(R-2));
      } else {
        printf("%5d %d %d %.10f %.10f ; %.10f %.10f bubble %.2f %d\n",
               iter, old, neww, max_dif, sum_dif/(L-2)/(R-2),
               v[neww][L/2][R/2], v[neww][L/3][R/3], bubble_volts, bvn);
      }
    }
    // check for convergence
    if (max_dif < 0.00000008) break;
 
  }
  printf(">> %d %.16f\n\n", iter, sum_dif);

  setup->bubble_volts = bubble_volts;
  setup->fully_depleted = 1;
  for (r=1; r<R; r++) {
    for (z=1; z<L; z++) {
      if (setup->point_type[z][r] < INSIDE) {
        undep[r][z] = ' ';
      } else if (undep[r][z] == 0) {
        undep[r][z] = '.';
      } else {
        if (undep[r][z] > 4) undep[r][z] = 'B';  // identifies pinch-off
        else undep[r][z] = '*';
        setup->fully_depleted = 0;
      }
    }
  }

  if (setup->vacuum_gap > 0) {   // restore impurity value along passivated surface
    for (r = 1; r < R; r++)
      setup->impurity[1][r] = setup->impurity[0][r];
  }

  return 0;
} /* ev_relax_undep */

/* -------------------------------------- wp_relax_undep ------------------- */
int wp_relax_undep(MJD_Siggen_Setup *setup) {
  int    old = 1, neww = 0, iter, r, z;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+2;
  int    R  = lrint(setup->xtal_radius/grid)+2;
  double eps_sum, v_sum, save_dif, pinched_sum1, pinched_sum2;
  double dif, sum_dif, max_dif;
  double ***v = setup->v, **eps_dr = setup->eps_dr, **eps_dz = setup->eps_dz;
  double ***dr = setup->dr, ***dz = setup->dz;
  double *s1 = setup->s1, *s2 = setup->s2;
  double e_over_E = 11.310; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0


  for (iter = 0; iter < setup->max_iterations; iter++) {

   double OR_fact = ((0.997 - 300.0/(L*R)) * (1.0 - 0.9/(double)(1+iter/6)));
    if (300.0/(L*R) > 0.5) OR_fact = (0.5 * (1.0 - 0.9/(double)(1+iter/6)));
    if (iter < 2) OR_fact = 0.0;

    old = neww;
    neww = 1 - neww;
    sum_dif = 0;
    max_dif = 0;
    pinched_sum1 = pinched_sum2 = 0.0;

    if (setup->vacuum_gap > 0) {   // modify impurity value along passivated surface
      for (r = 1; r < R; r++)      //   due to surface charge induced by capacitance
        setup->impurity[1][r] = -v[old][1][r] * 5.52e-4 * e_over_E * grid / setup->vacuum_gap;
    }

    /* start from z=1 and r=1 so that (z,r)=0 can be
       used for reflection symmetry around r=0 or z=0 */
    for (z = 1; z < L; z++) {
      /* manage r=0 reflection symmetry */
      setup->v[old][z][0] = setup->v[old][z][2];

      for (r = 1; r < R; r++) {
        if (setup->point_type[z][r] < INSIDE) continue;   // HV or point contact
        save_dif = v[old][z][r] - v[neww][z][r];      // step difference from previous iteration

        if (setup->point_type[z][r] < PINCHOFF) {       // normal bulk or passivated surface, no complications
          v_sum = (v[old][z+1][r] + v[old][z][r+1]*s1[r] +
                   v[old][z-1][r] + v[old][z][r-1]*s2[r]);
          if (r > 1) eps_sum = 4;
          else       eps_sum = 2 + s1[r] + s2[r];
        } else if (setup->point_type[z][r] == PINCHOFF) {  // in or adjacent to the ditch
          if (setup->point_type[z+1][r] < PINCHOFF) {
            pinched_sum1 += v[old][z+1][r]*eps_dz[z][r];
            pinched_sum2 += eps_dz[z][r];
          }
          if (setup->point_type[z][r+1] < PINCHOFF) {
            pinched_sum1 += v[old][z][r+1]*eps_dr[z][r]*s1[r];
            pinched_sum2 += eps_dr[z][r]*s1[r];
          }
          if (setup->point_type[z-1][r] < PINCHOFF) {
            pinched_sum1 += v[old][z-1][r]*eps_dz[z-1][r];
            pinched_sum2 += eps_dz[z-1][r];
          }
          if (setup->point_type[z][r-1] < PINCHOFF) {
            pinched_sum1 += v[old][z][r-1]*eps_dr[z][r-1]*s2[r];
            pinched_sum2 += eps_dr[z][r-1]*s2[r];
          }
          v_sum = pinched_sum1;
          eps_sum = pinched_sum2;
        } else if (setup->point_type[z][r] == CONTACT_EDGE) {  // adjacent to the contact
          v_sum = (v[old][z+1][r]*dz[1][z][r] + v[old][z][r+1]*dr[1][z][r] +
                   v[old][z-1][r]*dz[0][z][r] + v[old][z][r-1]*dr[0][z][r]);
          eps_sum = dz[1][z][r] + dr[1][z][r] + dz[0][z][r] + dr[0][z][r];
        } else if (setup->point_type[z][r] >= DITCH) {  // in or adjacent to the ditch
          v_sum = (v[old][z+1][r]*eps_dz[z  ][r] + v[old][z][r+1]*eps_dr[z][r  ]*s1[r] +
                   v[old][z-1][r]*eps_dz[z-1][r] + v[old][z][r-1]*eps_dr[z][r-1]*s2[r]);
          eps_sum = (eps_dz[z][r]   + eps_dr[z][r]  *s1[r] +
                     eps_dz[z-1][r] + eps_dr[z][r-1]*s2[r]);
        }

        if (setup->point_type[z][r] != PINCHOFF) {
          // calculate the interpolated mean potential and the effect of the space charge
          if (setup->vacuum_gap > 0 && z == 1)
            v[neww][z][r] = v_sum / eps_sum + setup->impurity[z][r];
          else
            v[neww][z][r] = v_sum / eps_sum;

          // calculate difference from last iteration, for convergence check
          dif = v[old][z][r] - v[neww][z][r];
          if (dif < 0) dif = -dif;
          sum_dif += dif;
          if (max_dif < dif) max_dif = dif;
          // do over-relaxation
          v[neww][z][r] += OR_fact*save_dif;
        }
      }
    }
    if (pinched_sum2 > 0.1) {
      for (z=1; z<L; z++) {
        for (r=1; r<R; r++) {
          if (setup->point_type[z][r] == PINCHOFF) {
            v[neww][z][r] = pinched_sum1 / pinched_sum2;
            dif = v[old][z][r] - v[neww][z][r];
            if (dif < 0) dif = -dif;
            sum_dif += dif;
            if (max_dif < dif) max_dif = dif;
          }
        }
      }
    }

    // report results for some iterations
    if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0) {
      printf("%5d %d %d %.10f %.10f ; %.10f %.10f\n",
             iter, old, neww, max_dif, sum_dif/(L-2)/(R-2),
             v[neww][L/2][R/2], v[neww][L/3][R/3]);
    }
    // check for convergence
    if (max_dif < 0.0000000001) break;

  }

  printf(">> %d %.16f\n\n", iter, sum_dif);

  return 0;
} /* wp_relax_undep */

/* -------------------------------------- interpolate ------------------- */
int interpolate(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup) {
  int    n, i, j, i2, j2, zmin, rmin, zmax, rmax;
  int    L  = lrint(old_setup->xtal_length/old_setup->xtal_grid)+3;
  int    R  = lrint(old_setup->xtal_radius/old_setup->xtal_grid)+3;
  int    L2 = lrint(setup->xtal_length/setup->xtal_grid)+3;
  int    R2 = lrint(setup->xtal_radius/setup->xtal_grid)+3;
  float  f, f1r, f1z, f2r, f2z;
  double ***v = setup->v, **ov = old_setup->v[1];

  /* the previous calculation was on a coarser grid...
     now copy/expand the potential to the new finer grid
  */
  n = (int) (old_setup->xtal_grid / setup->xtal_grid + 0.5);
  f = 1.0 / (float) n;
  printf("\ngrid %.4f -> %.4f; ratio = %d %.3f\n\n",
         old_setup->xtal_grid, setup->xtal_grid, n, f);
  for (i = i2 = 1; i < L-1; i++) {
    zmin = i2;
    zmax = i2 + n;
    if (zmax > L2-1) zmax = L2-1;
    for (j = j2 = 1; j < R-1; j++) {
      f1z = 0.0;
      rmin = j2;
      rmax = j2 + n;
      if (rmax > R2-1) rmax = R2-1;
      for (i2 = zmin; i2 < zmax; i2++) {
        f2z = 1.0 - f1z;
        f1r = 0.0;
        for (j2 = rmin; j2 < rmax; j2++) {
          f2r = 1.0 - f1r;
          v[0][i2][j2] = v[1][i2][j2] =      // linear interpolation
            f2z*f2r*ov[i][j  ] + f1z*f2r*ov[i+1][j  ] +
            f2z*f1r*ov[i][j+1] + f1z*f1r*ov[i+1][j+1];
          f1r += f;
        }
        f1z += f;
      }
      j2 = rmax;
    }
    i2 = zmax;
  }

  return 0;
} /* interpolate */


/* ************************* fields.h *************************** */
#define MAX_FNAME_LEN 512

static int nearest_field_grid_index(cyl_pt pt, cyl_int_pt *ipt, MJD_Siggen_Setup *setup);
static int grid_weights(cyl_pt pt, cyl_int_pt ipt, float out[2][2], MJD_Siggen_Setup *setup);
static cyl_pt efield(cyl_pt pt, cyl_int_pt ipt, MJD_Siggen_Setup *setup);
static int setup_efield(MJD_Siggen_Setup *setup);
static int setup_wp(MJD_Siggen_Setup *setup);
static int setup_velo(MJD_Siggen_Setup *setup);
static int efield_exists(cyl_pt pt, MJD_Siggen_Setup *setup);

/* field_setup
   given a field directory file, read electic field and weighting
   potential tables from files listed in directory
   returns 0 for success
*/
int field_setup(MJD_Siggen_Setup *setup){

  setup->rmin  = 0;
  setup->rmax  = setup->xtal_radius;
  setup->rstep = setup->xtal_grid;
  setup->zmin  = 0;
  setup->zmax  = setup->xtal_length;
  setup->zstep = setup->xtal_grid;
  if (setup->xtal_temp < MIN_TEMP) setup->xtal_temp = MIN_TEMP;
  if (setup->xtal_temp > MAX_TEMP) setup->xtal_temp = MAX_TEMP;

  TELL_NORMAL("rmin: %.2f rmax: %.2f, rstep: %.2f\n"
	      "zmin: %.2f zmax: %.2f, zstep: %.2f\n"
	      "Detector temperature is set to %.1f K\n",
	      setup->rmin, setup->rmax, setup->rstep,
	      setup->zmin, setup->zmax, setup->zstep,
	      setup->xtal_temp);

  if (setup_velo(setup) != 0){
    error("Failed to read drift velocity data from file: %s\n", 
	  setup->drift_name);
    return -1;
  }
  if (setup_efield(setup) != 0){
    error("Failed to read electric field data from file: %s\n", 
	  setup->field_name);
    return -1;
  }
  if (setup_wp(setup) != 0){
    error("Failed to read weighting potential from file %s\n",
	  setup->wp_name);
    return -1;
  }

  return 0;
}

static int efield_exists(cyl_pt pt, MJD_Siggen_Setup *setup){
  cyl_int_pt ipt;
  char ptstr[MAX_LINE];
  int  i, j, ir, iz;
  sprintf(ptstr, "(r,z) = (%.1f,%.1f)", pt.r, pt.z);
  if (outside_detector_cyl(pt, setup)){
    TELL_CHATTY("point %s is outside crystal\n", ptstr);
    return 0;
  }
  ipt.r = (pt.r - setup->rmin)/setup->rstep;  // CHECKED: no need for lrintf
  ipt.phi = 0;
  ipt.z = (pt.z - setup->zmin)/setup->zstep;  // CHECKED: no need for lrintf

  if (ipt.r < 0 || ipt.r + 1 >= setup->rlen ||
      ipt.z < 0 || ipt.z + 1 >= setup->zlen){
    TELL_CHATTY("point %s is outside wp table\n", ptstr);
    return 0;
  }
  for (i = 0; i < 2 ; i++){
    ir = ipt.r + i;
    for (j = 0; j < 2; j++){
      iz = ipt.z + j;
      if (setup->efld[ir][iz].r == 0.0 && setup->efld[ir][iz].z == 0.0) {
	TELL_CHATTY("point %s has no efield\n", ptstr);
	return 0;
      }
    }
  }
  TELL_CHATTY("point %s is in crystal\n", ptstr);
  return 1;
}

/* wpotential
   gives (interpolated) weighting potential at point pt, stored in wp
   returns 0 for success, 1 on failure
*/
int wpotential(point pt, float *wp, MJD_Siggen_Setup *setup){
  float w[2][2];
  int   i, j;
  cyl_int_pt ipt;
  cyl_pt cyl;

  // cyl = cart_to_cyl(pt);  // do not need to know phi, so save call to atan
  cyl.r = sqrt(pt.x*pt.x + pt.y*pt.y);
  cyl.z = pt.z;

  if (nearest_field_grid_index(cyl, &ipt, setup) < 0) return 1;
  grid_weights(cyl, ipt, w, setup);
  *wp = 0.0;
  for (i = 0; i < 2; i++){
    for (j = 0; j < 2; j++){
      *wp += w[i][j]*setup->wpot[ipt.r+i][ipt.z+j];
    }
  }

  return 0;
}

/* drift_velocity
   calculates drift velocity for charge q at point pt
   returns 0 on success, 1 on success but extrapolation was necessary,
   and -1 for failure
   anisotropic drift: crystal axes are assumed to be (x,y,z)
*/
int drift_velocity(point pt, float q, vector *velo, MJD_Siggen_Setup *setup){
  point  cart_en;
  cyl_pt e, en, cyl;
  cyl_int_pt ipt;
  int   i, sign;
  float abse, absv, f, a, b, c;
  float bp, cp, en4, en6;
  struct velocity_lookup *v_lookup1, *v_lookup2;

  /*  DCR: replaced this with faster code below, saves calls to atan and tan
  cyl = cart_to_cyl(pt);
  if (nearest_field_grid_index(cyl, &ipt, setup) < 0) return -1;
  e = efield(cyl, ipt, setup);
  abse = vector_norm_cyl(e, &en);
  en.phi = cyl.phi;
  cart_en = cyl_to_cart(en);
  */
  cyl.r = sqrt(pt.x*pt.x + pt.y*pt.y);
  cyl.z = pt.z;
  cyl.phi = 0;
  if (nearest_field_grid_index(cyl, &ipt, setup) < 0) return -1;
  e = efield(cyl, ipt, setup);
  abse = vector_norm_cyl(e, &en);
  if (cyl.r > 0.001) {
    cart_en.x = en.r * pt.x/cyl.r;
    cart_en.y = en.r * pt.y/cyl.r;
  } else {
    cart_en.x = cart_en.y = 0;
  }
  cart_en.z = en.z;

  /* find location in table to interpolate from */
  for (i = 0; i < setup->v_lookup_len - 2 && abse > setup->v_lookup[i+1].e; i++);
  v_lookup1 = setup->v_lookup + i;
  v_lookup2 = setup->v_lookup + i+1;
  f = (abse - v_lookup1->e)/(v_lookup2->e - v_lookup1->e);
  if (q > 0){
    a = (v_lookup2->ha - v_lookup1->ha)*f+v_lookup1->ha;
    b = (v_lookup2->hb- v_lookup1->hb)*f+v_lookup1->hb;
    c = (v_lookup2->hc - v_lookup1->hc)*f+v_lookup1->hc;
    bp = (v_lookup2->hbp- v_lookup1->hbp)*f+v_lookup1->hbp;
    cp = (v_lookup2->hcp - v_lookup1->hcp)*f+v_lookup1->hcp;
    setup->dv_dE = (v_lookup2->h100 - v_lookup1->h100)/(v_lookup2->e - v_lookup1->e);
  }else{
    a = (v_lookup2->ea - v_lookup1->ea)*f+v_lookup1->ea;
    b = (v_lookup2->eb- v_lookup1->eb)*f+v_lookup1->eb;
    c = (v_lookup2->ec - v_lookup1->ec)*f+v_lookup1->ec;
    bp = (v_lookup2->ebp- v_lookup1->ebp)*f+v_lookup1->ebp;
    cp = (v_lookup2->ecp - v_lookup1->ecp)*f+v_lookup1->ecp;
    setup->dv_dE = (v_lookup2->e100 - v_lookup1->e100)/(v_lookup2->e - v_lookup1->e);
  }
  /* velocity can vary from the direction of the el. field
     due to effect of crystal axes */
#define POW4(x) ((x)*(x)*(x)*(x))
#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
  en4 = POW4(cart_en.x) + POW4(cart_en.y) + POW4(cart_en.z);
  en6 = POW6(cart_en.x) + POW6(cart_en.y) + POW6(cart_en.z);
  absv = a + b*en4 + c*en6;
  sign = (q < 0 ? -1 : 1);
  setup->v_over_E = absv / abse;
  velo->x = sign*cart_en.x*(absv+bp*4*(cart_en.x*cart_en.x - en4)
			    + cp*6*(POW4(cart_en.x) - en6));
  velo->y = sign*cart_en.y*(absv+bp*4*(cart_en.y*cart_en.y - en4)
			    + cp*6*(POW4(cart_en.y) - en6));
  velo->z = sign*cart_en.z*(absv+bp*4*(cart_en.z*cart_en.z - en4)
			    + cp*6*(POW4(cart_en.z) - en6));
#undef POW4
#undef POW6
  return 0;
}

/* Find (interpolated or extrapolated) electric field for this point */
static cyl_pt efield(cyl_pt pt, cyl_int_pt ipt, MJD_Siggen_Setup *setup){
  cyl_pt e = {0,0,0}, ef;
  float  w[2][2];
  int    i, j;

  grid_weights(pt, ipt, w, setup);
  for (i = 0; i < 2; i++){
    for (j = 0; j < 2; j++){
      ef = setup->efld[ipt.r + i][ipt.z + j];
      e.r += ef.r*w[i][j];
      e.z += ef.z*w[i][j];
    }
  }
  e.phi = pt.phi;
  return e;
}


/* Find weights for 8 voxel corner points around pt for e/wp field*/
/* DCR: modified to work for both interpolation and extrapolation */
static int grid_weights(cyl_pt pt, cyl_int_pt ipt, float out[2][2],
			MJD_Siggen_Setup *setup){
  float r, z;

  r = (pt.r - setup->rmin)/setup->rstep - ipt.r;
  z = (pt.z - setup->zmin)/setup->zstep - ipt.z;

  out[0][0] = (1.0 - r) * (1.0 - z);
  out[0][1] = (1.0 - r) *        z;
  out[1][0] =        r  * (1.0 - z);
  out[1][1] =        r  *        z;
  return 0;
}


/*find existing integer field grid index closest to pt*/
/* added DCR */
static int nearest_field_grid_index(cyl_pt pt, cyl_int_pt *ipt,
				    MJD_Siggen_Setup *setup){
  /* returns <0 if outside crystal or too far from a valid grid point
              0 if interpolation is okay
              1 if we can find a point but extrapolation is needed
  */
  static cyl_pt  last_pt;
  static cyl_int_pt last_ipt;
  static int     last_ret = -99;
  cyl_pt new_pt;
  int    dr, dz;
  float  d[3] = {0.0, -1.0, 1.0};

  if (last_ret != -99 &&
      pt.r == last_pt.r && pt.z == last_pt.z) {
    *ipt = last_ipt;
    return last_ret;
  }
  last_pt = pt;
  last_ret = -2;

  if (outside_detector_cyl(pt, setup)) {
    last_ret = -1;
  } else{
    new_pt.phi = 0.0;
    for (dz=0; dz<3; dz++) {
      new_pt.z = pt.z + d[dz]*setup->zstep;
      for (dr=0; dr<3; dr++) {
	new_pt.r = pt.r + d[dr]*setup->rstep;
	if (efield_exists(new_pt, setup)) {
	  last_ipt.r = (new_pt.r - setup->rmin)/setup->rstep;  // CHECKED: do NOT use lrintf
	  last_ipt.phi = 0;
	  last_ipt.z = (new_pt.z - setup->zmin)/setup->zstep;  // CHECKED: do NOT use lrintf
	  *ipt = last_ipt;
	  if (dr == 0 && dz == 0) {
	    last_ret = 0;
	  } else {
	    last_ret = 1;
	  }
	  return last_ret;
	}
      }
    }
  }

  return last_ret;
}

/* setup_velo
   set up drift velocity calculations (read in table)
*/
static int setup_velo(MJD_Siggen_Setup *setup){
  static int vlook_sz = 0;
  static struct velocity_lookup *v_lookup;

  char  line[MAX_LINE], *c;
  FILE  *fp;
  int   i, v_lookup_len;
  struct velocity_lookup *tmp, v, v0;
  float sumb_e, sumc_e, sumb_h, sumc_h;

  double be=1.3e7, bh=1.2e7, thetae=200.0, thetah=200.0;  // parameters for temperature correction
  double pwre=-1.680, pwrh=-2.398, mue=5.66e7, muh=1.63e9; //     adopted for Ge   DCR Feb 2015
  double mu_0_1, mu_0_2, v_s_1, v_s_2, E_c_1, E_c_2, e, f;

  if (vlook_sz == 0) {
    vlook_sz = 10;
    if ((v_lookup = (struct velocity_lookup *)
	 malloc(vlook_sz*sizeof(*v_lookup))) == NULL) {
      error("malloc failed in setup_velo\n");
      return -1;
    }
  }
  if ((fp = fopen(setup->drift_name, "r")) == NULL){
    error("failed to open velocity lookup table file: '%s'\n", setup->drift_name);
    return -1;
  }
  line[0] = '#';
  c = line;
  while ((line[0] == '#' || line[0] == '\0') && c != NULL) c = fgets(line, MAX_LINE, fp);
  if (c == NULL) {
    error("Failed to read velocity lookup table from file: %s\n", setup->drift_name);
    fclose(fp);
    return -1;
  }
  TELL_CHATTY("Drift velocity table:\n"
	      "  e          e100    e110    e111    h100    h110    h111\n");   
  for (v_lookup_len = 0; ;v_lookup_len++){
    if (v_lookup_len == vlook_sz - 1){
      vlook_sz += 10;
      if ((tmp = (struct velocity_lookup *)
	   realloc(v_lookup, vlook_sz*sizeof(*v_lookup))) == NULL){
	error("realloc failed in setup_velo\n");
	fclose(fp);
	return -1;
      }
      v_lookup = tmp;
    }
    if (sscanf(line, "%f %f %f %f %f %f %f", 
	       &v_lookup[v_lookup_len].e,
	       &v_lookup[v_lookup_len].e100,
	       &v_lookup[v_lookup_len].e110,
	       &v_lookup[v_lookup_len].e111,
	       &v_lookup[v_lookup_len].h100,
	       &v_lookup[v_lookup_len].h110,
	       &v_lookup[v_lookup_len].h111) != 7){
      break; //assume EOF
    }	   
    //v_lookup[v_lookup_len].e *= 100; /*V/m*/
    tmp = &v_lookup[v_lookup_len];
    TELL_CHATTY("%10.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
		tmp->e, tmp->e100, tmp->e110, tmp->e111, tmp->h100, tmp->h110,tmp->h111);
    line[0] = '#';
    while ((line[0] == '#' || line[0] == '\0' ||
	    line[0] == '\n' || line[0] == '\r') && c != NULL) c = fgets(line, MAX_LINE, fp);
    if (c == NULL) break;
    if (line[0] == 'e' || line[0] == 'h') break; /* no more velocities data;
						    now reading temp correction data */
  }

  /* check for and decode temperature correction parameters */
  while (line[0] == 'e' || line[0] == 'h') {
    if (line[0] == 'e' &&
	sscanf(line+2, "%lf %lf %lf %lf", 
	       &mue, &pwre, &be, &thetae) != 4) break;//asume EOF
    if (line[0] == 'h' &&
	sscanf(line+2, "%lf %lf %lf %lf", 
	       &muh, &pwrh, &bh, &thetah) != 4) break;//asume EOF
    if (line[0] == 'e')
      TELL_CHATTY("electrons: mu_0 = %.2e x T^%.4f  B = %.2e  Theta = %.0f\n",
		  mue, pwre, be, thetae);
    if (line[0] == 'h')
      TELL_CHATTY("    holes: mu_0 = %.2e x T^%.4f  B = %.2e  Theta = %.0f\n",
		  muh, pwrh, bh, thetah);

    line[0] = '#';
    while ((line[0] == '#' || line[0] == '\0') && c != NULL) c = fgets(line, MAX_LINE, fp);
    if (c == NULL) break;
  }

  if (v_lookup_len == 0){
    error("Failed to read velocity lookup table from file: %s\n", setup->drift_name);
    return -1;
  }  
  v_lookup_len++;
  if (vlook_sz != v_lookup_len){
    if ((tmp = (struct velocity_lookup *) 
	 realloc(v_lookup, v_lookup_len*sizeof(*v_lookup))) == NULL){
      error("realloc failed in setup_velo. This should not happen\n");
      fclose(fp);
      return -1;
    }
    v_lookup = tmp;
    vlook_sz = v_lookup_len;
  }
  TELL_NORMAL("Drift velocity table has %d rows of data\n", v_lookup_len);
  fclose(fp);

  /*
    apply temperature dependence to mobilities;
    see drift_velocities.doc and tempdep.c
    The drift velocity reduces at higher temperature due to the increasing of
    scattering with the lattice vibration. We used a model by M. Ali Omar and
    L. Reggiani (Solid-State Electronics Vol. 30, No. 12 (1987) 1351) to
    calculate the temperature dependence.
  */
  /* electrons */
  TELL_NORMAL("Adjusting mobilities for temperature, from %.1f to %.1f\n", REF_TEMP, setup->xtal_temp);
  TELL_CHATTY("Index  field  vel_factor\n");
  mu_0_1 = mue * pow(REF_TEMP, pwre);
  v_s_1 = be * sqrt(tanh(0.5 * thetae / REF_TEMP));
  E_c_1 = v_s_1 / mu_0_1;
  mu_0_2 = mue * pow(setup->xtal_temp, pwre);
  v_s_2 = be * sqrt(tanh(0.5 * thetae / setup->xtal_temp));
  E_c_2 = v_s_2 / mu_0_2;
  for (i = 0; i < vlook_sz; i++){
    e = v_lookup[i].e;
    if (e < 1) continue;
    f = (v_s_2 * (e/E_c_2) / sqrt(1.0 + (e/E_c_2) * (e/E_c_2))) /
        (v_s_1 * (e/E_c_1) / sqrt(1.0 + (e/E_c_1) * (e/E_c_1)));
    v_lookup[i].e100 *= f;
    v_lookup[i].e110 *= f;
    v_lookup[i].e111 *= f;
    TELL_CHATTY("%2d %5.0f %f\n", i, e, f);
  }

  /* holes */
  mu_0_1 = muh * pow(REF_TEMP, pwrh);
  v_s_1 = bh * sqrt(tanh(0.5 * thetah / REF_TEMP));
  E_c_1 = v_s_1 / mu_0_1;
  mu_0_2 = muh * pow(setup->xtal_temp, pwrh);
  v_s_2 = bh * sqrt(tanh(0.5 * thetah / setup->xtal_temp));
  E_c_2 = v_s_2 / mu_0_2;
  for (i = 0; i < vlook_sz; i++){
    e = v_lookup[i].e;
    if (e < 1) continue;
    f = (v_s_2 * (e/E_c_2) / sqrt(1.0 + (e/E_c_2) * (e/E_c_2))) /
        (v_s_1 * (e/E_c_1) / sqrt(1.0 + (e/E_c_1) * (e/E_c_1)));
    v_lookup[i].h100 *= f;
    v_lookup[i].h110 *= f;
    v_lookup[i].h111 *= f;
    TELL_CHATTY("%2d %5.0f %f\n", i, e, f);
  }
  /* end of temperature correction */

  for (i = 0; i < vlook_sz; i++){
    v = v_lookup[i];
    v_lookup[i].ea =  0.5 * v.e100 -  4 * v.e110 +  4.5 * v.e111;
    v_lookup[i].eb = -2.5 * v.e100 + 16 * v.e110 - 13.5 * v.e111;
    v_lookup[i].ec =  3.0 * v.e100 - 12 * v.e110 +  9.0 * v.e111;
    v_lookup[i].ha =  0.5 * v.h100 -  4 * v.h110 +  4.5 * v.h111;
    v_lookup[i].hb = -2.5 * v.h100 + 16 * v.h110 - 13.5 * v.h111;
    v_lookup[i].hc =  3.0 * v.h100 - 12 * v.h110 +  9.0 * v.h111;
  }
  v_lookup[0].ebp = v_lookup[0].ecp = v_lookup[0].hbp = v_lookup[0].hcp = 0.0;
  sumb_e = sumc_e = sumb_h = sumc_h = 0.0;
  for (i = 1; i < vlook_sz; i++){
    v0 = v_lookup[i-1];
    v = v_lookup[i];
    sumb_e += (v.e - v0.e)*(v0.eb+v.eb)/2;
    sumc_e += (v.e - v0.e)*(v0.ec+v.ec)/2;
    sumb_h += (v.e - v0.e)*(v0.hb+v.hb)/2;
    sumc_h += (v.e - v0.e)*(v0.hc+v.hc)/2;
    v_lookup[i].ebp = sumb_e/v.e;
    v_lookup[i].ecp = sumc_e/v.e;
    v_lookup[i].hbp = sumb_h/v.e;
    v_lookup[i].hcp = sumc_h/v.e;
  }

  setup->v_lookup = v_lookup;
  setup->v_lookup_len = v_lookup_len;

  return 0;
}


/* This may or may not break if we switch to a non-integer grid*/
/*setup_efield
  read electric field data from file, apply sanity checks
  returns 0 for success
*/
static int setup_efield(MJD_Siggen_Setup *setup){
  FILE   *fp;
  char   line[MAX_LINE], *cp;
  int    i, j, lineno = 0;
  float  v, eabs, er, ez;
  cyl_pt cyl, **efld;

  if ((fp = fopen(setup->field_name, "r")) == NULL){
    error("failed to open electric field table: %s\n", setup->field_name);
    return 1;
  }
  
  setup->rlen = lrintf((setup->rmax - setup->rmin)/setup->rstep) + 1;
  setup->zlen = lrintf((setup->zmax - setup->zmin)/setup->zstep) + 1;
  TELL_CHATTY("rlen, zlen: %d, %d\n", setup->rlen, setup->zlen);

  // here I assume that r, zlen never change from their initial values, which is reasonable
  if ((efld = (cyl_pt **) malloc(setup->rlen*sizeof(*efld))) == NULL) {
    error("Malloc failed in setup_efield\n");
    fclose(fp);
    return 1;
  }
  for (i = 0; i < setup->rlen; i++){
    if ((efld[i] = (cyl_pt *) malloc(setup->zlen*sizeof(*efld[i]))) == NULL){
      error("Malloc failed in setup_efield\n");
      //NB: potential memory leak here.
      fclose(fp);
      return 1;
    }
    memset(efld[i], 0, setup->zlen*sizeof(*efld[i]));
  }
  TELL_NORMAL("Reading electric field data from file: %s\n", setup->field_name);

  /*now read the table*/
  while(fgets(line, MAX_LINE, fp) != NULL){
    lineno++;
    for (cp = line; isspace(*cp) && *cp != '\0'; cp++);
    if (*cp == '#' || !strlen(cp)) continue;
    if (sscanf(line, "%f %f %f %f %f %f", 
	       &cyl.r, &cyl.z, &v, &eabs, &er, &ez) != 6){
      error("failed to read electric field data from line no %d\n"
	    "of file %s\n", lineno, setup->field_name);
      fclose(fp);
      return 1;
    }
    i = lrintf((cyl.r - setup->rmin)/setup->rstep);
    j = lrintf((cyl.z - setup->zmin)/setup->zstep);
    if (i < 0 || i >= setup->rlen || j < 0 || j >= setup->zlen) {
      // error("Error in efield line %d, i = %d, j = %d\n", line, i, j);
      continue;
    }
    cyl.phi = 0;
    if (outside_detector_cyl(cyl, setup)) continue;
    efld[i][j].r = er;
    efld[i][j].z = ez;
    efld[i][j].phi = 0;
  }      

  TELL_NORMAL("Done reading %d lines of electric field data\n", lineno);
  fclose(fp);

  setup->efld = efld;
  for (i = 0; i < setup->rlen; i++) setup->efld[i] = efld[i];

  return 0;
}

/*setup_wp
  read weighting potential values from files. returns 0 on success*/
static int setup_wp(MJD_Siggen_Setup *setup){
  FILE   *fp;
  char   line[MAX_LINE], *cp;
  int    i, j, lineno;
  cyl_pt cyl;
  float  wp, **wpot;

  setup->rlen = lrintf((setup->rmax - setup->rmin)/setup->rstep) + 1;
  setup->zlen = lrintf((setup->zmax - setup->zmin)/setup->zstep) + 1;
  TELL_CHATTY("rlen, zlen: %d, %d\n", setup->rlen, setup->zlen);

  //assuming rlen, zlen never change as for setup_efld
  if ((wpot = (float **) malloc(setup->rlen*sizeof(*wpot))) == NULL){
    error("Malloc failed in setup_wp\n");
    return 1;
  }
  for (i = 0; i < setup->rlen; i++){
    if ((wpot[i] = (float *) malloc(setup->zlen*sizeof(*wpot[i]))) == NULL){  
      error("Malloc failed in setup_wp\n");
      //NB: memory leak here.
      return 1;
    }
    memset(wpot[i], 0, setup->zlen*sizeof(*wpot[i]));
  }
  if ((fp = fopen(setup->wp_name, "r")) == NULL){
    error("failed to open file: %s\n", setup->wp_name);
    return -1;
  }
  lineno = 0;
  TELL_NORMAL("Reading weighting potential from file: %s\n", setup->wp_name);
  while (fgets(line, MAX_LINE, fp) != NULL){
    lineno++;
    for (cp = line; isspace(*cp) && *cp != '\0'; cp++);
    if (*cp == '#' || !strlen(cp)) continue;
    if (sscanf(line, "%f %f %f\n",&cyl.r, &cyl.z, &wp) != 3){ 
      error("failed to read weighting potential from line %d\n"
	    "line: %s", lineno, line);
      fclose(fp);
      return 1;
    }
    i = lrintf((cyl.r - setup->rmin)/setup->rstep);
    j = lrintf((cyl.z - setup->zmin)/setup->zstep);
    if (i < 0 || i >= setup->rlen || j < 0 || j >= setup->zlen) continue;
    if (outside_detector_cyl(cyl, setup)) continue;
    wpot[i][j] = wp;
  }
  TELL_NORMAL("Done reading %d lines of WP data\n", lineno);
  fclose(fp);

  setup->wpot = wpot;
  for (i = 0; i < setup->rlen; i++) setup->wpot[i] = wpot[i];

  return 0;
}


/* free malloc()'ed memory and do other cleanup*/
int fields_finalize(MJD_Siggen_Setup *setup){
  int i;

  for (i = 0; i < lrintf((setup->rmax - setup->rmin)/setup->rstep) + 1; i++){
    free(setup->efld[i]);
    free(setup->wpot[i]);
  }
  free(setup->efld);
  free(setup->wpot);
  free(setup->v_lookup);
  setup->efld = NULL;
  setup->wpot = NULL;
  setup->v_lookup = NULL;

  return 1;
}

void set_temp(float temp, MJD_Siggen_Setup *setup){
  if (temp < MIN_TEMP || temp > MAX_TEMP){
    error("temperature out of range: %f\n", temp);
  }else{
    setup->xtal_temp = temp;
    error("temperature set to %f\n", temp);
    /* re-read velocities and correct them to the new temperature value */
    setup_velo(setup);
  }
}

/* ************************* calc_signal.h *************************** */
#define HOLE_CHARGE 1.0
#define ELECTRON_CHARGE -1.0
/* the following is the diffusion coefficient for holes in Ge at 77K
               at low field (~ 100 V/cm) 
   the diffusion coefficient drops at higher fields, and higher temperatures
   see Jacoboni et al., Phys. Rev. B24, 2 (1981) 1014-1026.
   size sigma = sqrt(2Dt), t = time, D = mu*k*T/e
   mu = mobility, k = Boltzmann const., T = temp, e = electron charge
   mu_h = 4e4 cm^2/V/s, mu_e = 5e4 cm^2/V/s at 77K, so
   D_h = 265 cm^2/s, D_e = 332 cm^2/s
   and goes down roughly as 1/Temp (since mu goes as T^-1.7 or T^-2.3)

   we also convert (2Dt) from sigma-squared to FWHM-squared

   for Si at 300K, 
   mu_h = 450 cm^2/V/s, mu_e = 1500 cm^2/V/s, so
   D_h = 12 cm^2/s, D_e = 39 cm^2/s
*/
/*  here are some definitions used for an old method, where I calculated FWHM-squared:
// germanium:  2Dt = sigma^2; we want  FWHM^2 in mm^2 / ns
#define TWO_TIMES_DIFFUSION_COEF_H \
        (2.0 * 2.355*2.355 * 2.65e-5 * setup->step_time_calc * 77.0/setup->xtal_temp)
#define TWO_TIMES_DIFFUSION_COEF_E \
        (2.0 * 2.355*2.355 * 3.32e-5 * setup->step_time_calc * 77.0/setup->xtal_temp)
// silicon:
#define TWO_TIMES_DIFFUSION_COEF_H_Si \
        (1.3e-5 * setup->step_time_calc * 300.0/setup->xtal_temp)
#define TWO_TIMES_DIFFUSION_COEF_E_Si \
        (4.3e-5 * setup->step_time_calc * 300.0/setup->xtal_temp)
*/
/* In the new method, I use dsigma/dt = D/sigma to calculate FWHM */
#define DIFFUSION_COEF   (setup->v_over_E * 0.67)
/* above is my own approximate parameterization of measurements of Jacoboni et al.
   0.67 = 2.355 * 2.355 * 0.12    to get D in mm2/s, and scaled to FWHM2/sigma2
   v_over_E = drift velocity / electric field   ~  mu
   note that Einstein's equation is D = mu*kT/e
   kT/e ~ 0.007/V ~ 0.07 mm/Vcm, => close enough to 0.12, okay
 */

/* prototypes for module-private functions*/
//static int make_signal(point pt, float *signal, float q, MJD_Siggen_Setup *setup);
//static float charge_trapping(vector dx, float q); //trapping

/* signal_calc_init
   read setup from configuration file,
   then read the electric field and weighting potential,
   and initialize the signal calculation variables
   returns 0 for success
*/
int signal_calc_init(char *config_file_name, MJD_Siggen_Setup *setup) {

  if (read_config(config_file_name, setup)) return 1;

  TELL_CHATTY("r: %.2f  z: %.2f\n", setup->xtal_radius, setup->xtal_length);
  setup->ntsteps_out = setup->time_steps_calc /
    lrintf(setup->step_time_out/setup->step_time_calc);
  TELL_NORMAL("Will use %d time steps in calculations, each %.2f ns long;\n"
	      "the output signals will have %d time steps, each %.2f ns long\n", 
	      setup->time_steps_calc, setup->step_time_calc,
	      setup->ntsteps_out, setup->step_time_out);

  TELL_NORMAL("Reading field data...\n");
  if (field_setup(setup) != 0) return -1;
  
  if ((setup->dpath_e = (point *) malloc(setup->time_steps_calc*sizeof(point))) == NULL ||
      (setup->dpath_h = (point *) malloc(setup->time_steps_calc*sizeof(point))) == NULL) {
    error("Path malloc failed\n");
    return -1;
  }

  tell("Setup of signal calculation done\n");
  return 0;
}

/* get_signal
   calculate signal for point pt. Result is placed in signal_out array
   returns -1 if outside crystal
   if signal_out == NULL => no signal is stored
*/
int get_signal(point pt, float *signal_out, MJD_Siggen_Setup *setup) {
  static float *signal, *sum, *tmp;
  static int tsteps = 0;
  float w, x, y;
  char  tmpstr[MAX_LINE];
  int   j, k, l, dt, err, comp_f;

  /* first time -- allocate signal and sum arrays */
  if (tsteps != setup->time_steps_calc) {
    tsteps = setup->time_steps_calc;
    if ((signal = (float *) malloc(tsteps*sizeof(*signal))) == NULL ||
        (tmp    = (float *) malloc(tsteps*sizeof(*tmp))) == NULL ||
        (sum    = (float *) malloc(tsteps*sizeof(*sum))) == NULL) {
      error("malloc failed in get_signal\n");
      return -1;
    }
  }

  for (j = 0; j < tsteps; j++) signal[j] = 0.0;

  if (outside_detector(pt, setup)) {
    TELL_CHATTY("Point %s is outside detector!\n", pt_to_str(tmpstr, MAX_LINE, pt));
    return -1;
  }
  TELL_CHATTY("Calculating signal for %s...\n", pt_to_str(tmpstr, MAX_LINE, pt));

  memset(setup->dpath_e, 0, tsteps*sizeof(point));
  memset(setup->dpath_h, 0, tsteps*sizeof(point));

  err = make_signal(pt, signal, ELECTRON_CHARGE, setup);
  err = make_signal(pt, signal, HOLE_CHARGE, setup);
  /* make_signal returns 0 for success; require hole signal but not electron */

  /* change from current signal to charge signal, i.e.
     each time step contains the summed signals of all previous time steps */
  for (j = 1; j < tsteps; j++) signal[j] += signal[j-1];

  if (signal_out != NULL) {

    if (setup->charge_cloud_size > 0.001 || setup->use_diffusion) {
      /* convolute with a Gaussian to correct for charge cloud size
	 and initial velocity
	 charge_cloud_size = initial FWHM of charge cloud, in mm,
	 NOTE this uses initial velocity of holes only;
	 this may not be quite right if electron signal is strong */
      /* difference in time between center and edge of charge cloud */
      dt = (int) (1.5f + setup->charge_cloud_size /
		  (setup->step_time_calc * setup->initial_vel));
      if (setup->initial_vel < 0.00001f) dt = 0;
      TELL_CHATTY("Initial vel, size, dt = %f mm/ns, %f mm, %d steps\n",
		  setup->initial_vel, setup->charge_cloud_size, dt);
      if (setup->use_diffusion) {
	dt = (int) (1.5f + setup->final_charge_size /
		    (setup->step_time_calc * setup->final_vel));
	TELL_CHATTY("  Final vel, size, dt = %f mm/ns, %f mm, %d steps\n",
		    setup->final_vel, setup->final_charge_size, dt);
      }
      if (dt > 1) {
	/* Gaussian */
	w = ((float) dt) / 2.355;
	l = dt/10;     // use l to speed up convolution of waveform with gaussian;
	if (l < 1) {   // instead of using every 1-ns step, use steps of FWHM/10
	  l = 1;
	} else if (setup->step_time_out > setup->preamp_tau) {
	  if (l > setup->step_time_out/setup->step_time_calc)
	    l = setup->step_time_out/setup->step_time_calc;
	} else {
	  if (l > setup->preamp_tau/setup->step_time_calc)
	    l = setup->preamp_tau/setup->step_time_calc;
	}
	// TELL_CHATTY(">> l: %d\n", l);
	for (j = 0; j < tsteps; j++) {
	  sum[j] = 1.0;
	  tmp[j] = signal[j];
	}
	for (k = l; k < 2*dt; k+=l) {
	  x = ((float) k)/w;
	  y = exp(-x*x/2.0);
	  for (j = 0; j < tsteps - k; j++){
	    sum[j] += y;
	    tmp[j] += signal[j+k] * y;
	    sum[j+k] += y;
	    tmp[j+k] += signal[j] * y;
	  }
	}
        for (j = 0; j < tsteps; j++){
          signal[j] = tmp[j]/sum[j];
        }
      }
    }

    /* now, compress the signal and place it in the signal_out array;
       truncate the signal if time_steps_calc % ntsteps_out != 0 */
    comp_f = setup->time_steps_calc/setup->ntsteps_out;
    for (j = 0; j < setup->ntsteps_out; j++) signal_out[j] = 0;
    for (j = 0; j < setup->ntsteps_out*comp_f; j++)
      signal_out[j/comp_f] += signal[j]/comp_f;

    /* do RC integration for preamp risetime */
    if (setup->preamp_tau/setup->step_time_out >= 0.1f)
      rc_integrate(signal_out, signal_out,
		   setup->preamp_tau/setup->step_time_out, setup->ntsteps_out);
  }

  /* make_signal returns 0 for success; require hole signal but not electron */
  if (err) return -1;
  return 1;
}

/* make_signal
   Generates the signal originating at point pt, for charge q
   returns 0 for success
*/
int make_signal(point pt, float *signal, float q, MJD_Siggen_Setup *setup) {
  float  wpot, wpot2=0, dwpot;
  char   tmpstr[MAX_LINE];
  point  new_pt;
  vector v, dx;
  float  vel0, vel1 = 0, wpot_old=-1;
  // double diffusion_coeff;
  double repulsion_fact = 0.0, ds2, ds3, dv, ds_dt;
  int    ntsteps, i, t, n, collect2pc, low_field=0, surface_drift=0, stop_drifting = 0;

  new_pt = pt;
  collect2pc = ((q > 0 && setup->impurity_z0 < 0) ||  // holes for p-type 
		(q < 0 && setup->impurity_z0 > 0));   // electrons for n-type
  /*
  if (q > 0) {
    diffusion_coeff = TWO_TIMES_DIFFUSION_COEF_H;
  } else {
    diffusion_coeff = TWO_TIMES_DIFFUSION_COEF_E;
  }
  */
  ntsteps = setup->time_steps_calc;
  for (t = 0; drift_velocity(new_pt, q, &v, setup) >= 0 && !stop_drifting; t++) { 
    if (q > 0) {
      setup->dpath_h[t] = new_pt;
    } else {
      setup->dpath_e[t] = new_pt;
    }
    if (collect2pc) {
      if (t == 0) {
	vel1 = setup->final_vel = setup->initial_vel = vector_length(v);
	setup->final_charge_size = setup->charge_cloud_size;
	if (setup->use_diffusion) {
	  if (setup->final_charge_size < 0.01) setup->final_charge_size = 0.01;
	  /* for a spherically symmetric charge cloud, the equivalent
	     delta-E at a distance of 1 sigma from the cloud center is
	     dE = Q/(4*pi*epsilon*sigma^2)  (Q is charge inside the 3D 1-sigma envelope)
	     dE (V/cm) = Q (C) * 1/(4*pi*epsilon) (N m2 / C2) / sigma2 (mm2)
	     1 V/m = 1 N/C
	     dE (V/cm) = Q (C) * 1/(4*pi*epsilon) (V m / C) / sigma2 (mm2)
	     dE (V/cm) = repulsion_fact * FWHM/sigma / (FWHM^2) (mm2), so
	     repulsion_fact = (FWHM/sigma)^3 * Q (C) * 1/(4*pi*epsilon) (V m / C) * mm/m * mm/cm
	  */
	  if (setup->energy > 0.1) {  // set up charge cloud self-repulsion
	    repulsion_fact = setup->energy * 0.67*0.67*0.67 / 0.003; // charge in 1 sigma (3D)
	    repulsion_fact /= 6.241e18;        // convert to Coulombs
	    repulsion_fact *= 9.0e13/16.0;     // 1/(4*pi*epsilon)  (N m2 / C2) * 1e4
	    repulsion_fact *= 2.355*2.355*2.355;      // convert FWHM to sigma
	  }
	}
	TELL_CHATTY("initial v: %f (%e %e %e)\n",
		    setup->initial_vel, v.x, v.y, v.z);
      } else if (setup->use_diffusion) {
	vel0 = vel1;
	vel1 = vector_length(v);
	setup->final_charge_size *= vel1/vel0;  // effect of acceleration
	// include effects of acceleration and diffusion on cloud size
	dv = repulsion_fact * setup->dv_dE /        // effect of repulsion
	        (setup->final_charge_size*setup->final_charge_size);
	// FIXME? this next line could more more fine-grained
	if (dv > 0.05) dv = 0.05;  // on account of drift velocity saturation
	ds_dt = dv + DIFFUSION_COEF/setup->final_charge_size;  // effect of diffusion
	if (ds_dt > 0.05 || ds_dt * setup->step_time_calc > 0.1) {
	  // nonlinear growth due to small size; need more careful calculation
	  TELL_CHATTY("ds_dt = %.2f; size = %.2f", ds_dt, setup->final_charge_size);
	  // ds_dt = 0.05;  // artificially limit nonlinear growth
	  ds2 = 2.0 * DIFFUSION_COEF * setup->step_time_calc; // increase^2 from diff.
	  ds3 = (setup->final_charge_size*setup->final_charge_size *
		 (setup->final_charge_size +
		  3.0 * dv * setup->step_time_calc));         // FWHM^3 after repulsion
	  setup->final_charge_size = sqrt(ds2 + pow(ds3, 0.6667)); 
	  TELL_CHATTY(" -> %.2f\n", setup->final_charge_size);
	} else {
	  setup->final_charge_size +=  ds_dt * setup->step_time_calc;  // effect of diff. + rep.
	}
      }
    }

    TELL_CHATTY("pt: (%.2f %.2f %.2f), v: (%e %e %e)",
		new_pt.x, new_pt.y, new_pt.z, v.x, v.y, v.z);
    if (0 && t >= ntsteps - 2) {   // DRC removed (if(0)) Oct 2019; t>ntsteps now dealt with below
      if (collect2pc || wpot > WP_THRESH_ELECTRONS) {
	/* for p-type, this is hole or electron+high wp */
	TELL_CHATTY("\nExceeded maximum number of time steps (%d)\n", ntsteps);
	low_field = 1;
	// return -1;
      }
      break;
    }
    if (wpotential(new_pt, &wpot, setup) != 0) {
      TELL_NORMAL("\nCan calculate velocity but not WP at %s!\n",
		  pt_to_str(tmpstr, MAX_LINE, new_pt));
      return -1;
    }
    if (wpot < 0.0) wpot = 0.0;
    TELL_CHATTY(" -> wp: %.4f\n", wpot);

    /* ------------- DCR added Oct 2019: if WP is very small or large, then stop drifting */
    if (!collect2pc &&    wpot < 5.0e-5) stop_drifting = 2;    // drifting to outside
    if (collect2pc && 1.0-wpot < 5.0e-5) stop_drifting = 3;    // drifting to point contact
    if (t >= setup->time_steps_calc - 2) stop_drifting = 1;    // have run out of time...

    if (t > 0) signal[t] += q*(wpot - wpot_old);
    // FIXME! Hack added by DCR to deal with undepleted point contact
    if (wpot >= 0.999 && (wpot - wpot_old) < 0.0002) {
      low_field = 1;
      break;
    }
    wpot_old = wpot;

    dx = vector_scale(v, setup->step_time_calc);
    if (surface_drift && dx.z < 0) {
      dx.x *= setup->surface_drift_vel_factor;  // Hmmm... should the default be zero or one?
      dx.y *= setup->surface_drift_vel_factor;
      dx.z = 0;
    }
    new_pt = vector_add(new_pt, dx);
    // q = charge_trapping(dx, q); //FIXME

    // look for charges on passivated surface of a PPC detector
    if (new_pt.z < 0 &&                                    // at or below surface, and
        (setup->wrap_around_radius <= setup->pc_radius ||  // this is a PPC detector
         new_pt.x*new_pt.x + new_pt.y*new_pt.y <           // or point is inside wrap-around
         setup->wrap_around_radius*setup->wrap_around_radius)) {
      TELL_CHATTY(" -> Passivated surface! q = %.2f  r = %.2f\n",
                  q, sqrt(new_pt.x*new_pt.x + new_pt.y*new_pt.y));
      //break;
      surface_drift = 1;
      new_pt.z = 0;
    }

  }
  if (t == 0) {
    TELL_CHATTY("The starting point %s is outside the WP or field.\n",
		pt_to_str(tmpstr, MAX_LINE, pt));
    return -1;
  }

  if (low_field) {
    TELL_CHATTY("Low field near point contact; this may or may not be a problem.\n");
  } else {
    TELL_CHATTY("Drifted to edge of WP or field grid, point: %s q: %.2f\n", 
		pt_to_str(tmpstr, MAX_LINE, new_pt), q);
  }
  if (!low_field && stop_drifting<2) {
    /* figure out how much we must drift to get to the crystal boundary */
    for (n = 0; n+t < ntsteps && !outside_detector(new_pt, setup); n++){
      new_pt = vector_add(new_pt, dx);
      if (q > 0) setup->dpath_h[t+n] = new_pt;
      else setup->dpath_e[t+n] = new_pt;
    }
    if (n == 0) n = 1; /* always drift at least one more step */
    // TELL_CHATTY(
    TELL_NORMAL("q: %.1f t: %d n: %d ((%.2f %.2f %.2f)=>(%.2f %.2f %.2f))\n", 
		q, t, n, pt.x, pt.y, pt.z, new_pt.x, new_pt.y, new_pt.z);

    if (n + t >= ntsteps){
      n = ntsteps - t;
      if (q > 0 || wpot > WP_THRESH_ELECTRONS) { /* hole or electron+high wp */
	TELL_CHATTY("Exceeded maximum number of time steps (%d)\n", ntsteps);
        /* check WP to see if we have produced most of the signal */
        if ((wpot < 0.95 || wpot > 0.05) &&
            wpotential(new_pt, &wpot2, setup) != 0) {
          TELL_CHATTY("Cannot finish drifting to make at least 95% of signal.\n");
          return -1;  /* FIXME: could this be improved? */
        }
        /* drift to new_pt and wpot2 */
        dwpot = (wpot2 - wpot)/n;
      }
    } else {
      /* make WP go gradually to 1 or 0 */
      if (wpot > 0.3) {
        dwpot = (1.0 - wpot)/n;
      } else {
        dwpot = - wpot/n;
      }
    }

    /* now drift the final n steps */
    dx = vector_scale(v, setup->step_time_calc);
    if (new_pt.z > 0) {               // charges NOT on passivated surface
      for (i = 0; i < n; i++){
        signal[i+t] += q*dwpot;
        // q = charge_trapping(dx, q); //FIXME
      }
    }
  }
  TELL_CHATTY("q:%.2f pt: %s\n", q, pt_to_str(tmpstr, MAX_LINE, pt));
  if (q > 0) setup->final_vel = vector_length(v);

  return 0;
}

//FIXME -- placeholder function. Even parameter list is dodgy
/*
static float charge_trapping(vector dx, float q){
  return q;
}
*/

int rc_integrate(float *s_in, float *s_out, float tau, int time_steps){
  int   j;
  float s_in_old, s;  /* DCR: added so that it's okay to
			 call this function with s_out == s_in */
  
  if (tau < 1.0f) {
    for (j = time_steps-1; j > 0; j--) s_out[j] = s_in[j-1];
    s_out[0] = 0.0;
  } else {
    s_in_old = s_in[0];
    s_out[0] = 0.0;
    for (j = 1; j < time_steps; j++) {
      s = s_out[j-1] + (s_in_old - s_out[j-1])/tau;
      s_in_old = s_in[j];
      s_out[j] = s;
    }
  }
  return 0;
}

/* signal_calc_finalize
 * Clean up (free arrays, close open files...)
 */
int signal_calc_finalize(MJD_Siggen_Setup *setup){
  fields_finalize(setup);
  free(setup->dpath_h);
  free(setup->dpath_e);
  return 0;
}

int drift_path_e(point **pp, MJD_Siggen_Setup *setup){
  *pp = setup->dpath_e;
  return setup->time_steps_calc;
}
int drift_path_h(point **pp, MJD_Siggen_Setup *setup){
  *pp = setup->dpath_h;
  return setup->time_steps_calc;
}

/* tell
   write to stdout, provided that verb_level is above the threshold */
void tell(const char *format, ...){
  va_list ap;

  va_start(ap, format);
  vprintf(format, ap);
  va_end(ap);
  return;
}

/*error
  report error messages to stderr */
void error(const char *format, ...) {
  va_list ap;

  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  return;
}


namespace py = pybind11;
/* ************************* python interfaces *************************** */
py::list run_get_signal( point pt, MJD_Siggen_Setup *setup){
  // float *signal_out;
  // MJD_Siggen_Setup setup;
  py::list res;
  static float *signal_out;

  // signal_calc_init(config_file_name, &setup);

  if (signal_out == NULL){//first call
    if ((signal_out = (float *) malloc(setup->ntsteps_out*sizeof(*signal_out))) == NULL){
      printf("Malloc failed\n");
      return res;
    }
  }
  get_signal(pt, signal_out, setup);
  printf("asd3\n");
  // for (auto signal_out : li){
  //   res::append(li);
  // }
  // int arrSize = sizeof(signal_out)/sizeof(signal_out[0]);
  // printf("arrSize %d, %d", arrSize, setup.ntsteps_out);
  for( size_t i = 0 ; i < setup->ntsteps_out ; i++ )
  {
    res.append(signal_out[i]);
  }
  return res;
}



/* ************************* python interfaces *************************** */



// _test имя нашего модуля
PYBIND11_MODULE(mjd_fieldgen, m) {
    // m.def("main", &main);

    m.def("run_mjd_fieldgen", &main, "Run mjd_fieldgen",
      py::arg("config_file") = "PPC.config",
      py::arg("w") = -1, 
      py::arg("p") = -1, 
      py::arg("d") = -1,
      py::arg("b") = 0, 
      py::arg("impurity_profile") = "");


    py::class_<point>(m, "point")
        .def(py::init())
        .def_readwrite("x", &point::x)
        .def_readwrite("y", &point::y)
        .def_readwrite("z", &point::z)
        .def("__repr__",
          [](const point &a) {
              return "(x: "+ std::to_string(a.x)+", y: "+std::to_string(a.y)+", z: "+std::to_string(a.z) + ")";
              // return "asd";
        });
    py::class_<cyl_pt>(m, "cyl_pt")
        .def(py::init())
        .def_readwrite("r", &cyl_pt::r)
        .def_readwrite("phi", &cyl_pt::phi)
        .def_readwrite("z", &cyl_pt::z)
        .def("__repr__",
          [](const cyl_pt &a) {
              return "(r: "+ std::to_string(a.r)+", phi: "+std::to_string(a.phi)+", z: "+std::to_string(a.z) + ")";
              // return "asd";
        });
    py::class_<MJD_Siggen_Setup>(m, "MJD_Siggen_Setup")
        .def(py::init())
        .def_readwrite("time_steps_calc", &MJD_Siggen_Setup::time_steps_calc)
        .def_readwrite("xtal_length", &MJD_Siggen_Setup::xtal_length)
        .def_readwrite("xtal_radius", &MJD_Siggen_Setup::xtal_radius);
    m.def("run_get_signal", &run_get_signal, "Run asdf",
      // py::arg("config_file") = "PPC.config",
      py::arg("pt"),
      py::arg("setup"));
    m.def("signal_calc_init", &signal_calc_init, "asdfa",
      py::arg("config_file_name") = "PPC.config",
      py::arg("setup"));
    m.def("outside_detector_cyl", &outside_detector_cyl, "rasdfa",
      py::arg("pt"),
      py::arg("setup"));
    m.def("efield_exists", &efield_exists, "rasdfa",
      py::arg("pt"),
      py::arg("setup"));
};