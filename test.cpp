#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
#include <cmath>
#include <elina_abstract0.h>
#include <opt_oct.h>

static int shouldClear = 1;

double norm(std::vector<double> &v1, std::vector<double> &v2) {
  std::vector<double> v3;
  for (int i = 0; i < v1.size(); i++) {
    v3.push_back(v2[i] - v1[i]);
  }
  return std::sqrt(std::inner_product(v3.begin(), v3.end(), v3.begin(), 0));
}

std::vector<double> compute_center(elina_manager_t* man, elina_abstract0_t* abs) {
  elina_interval_t** itv = elina_abstract0_to_box(man, abs);
  int dims = elina_abstract0_dimension(man, abs).realdim;
  std::vector<double> center(dims);
  for (int i = 0; i < dims; i++) {
    double l, u;
    elina_double_set_scalar(&l, itv[i]->inf, MPFR_RNDN);
    elina_double_set_scalar(&u, itv[i]->sup, MPFR_RNDN);
    center[i] = (l + u) / 2.0;
  }

  for (int i = 0; i < dims; i++) {
    elina_interval_free(itv[i]);
  }
  free(itv);

  return center;
}


class Powerset {
  public:
    int size;
    elina_manager_t *man;
    std::vector<elina_abstract0_t*> disjuncts;
    std::vector<std::vector<double>> centers;
    Powerset(const Powerset &);
    Powerset(std::vector<elina_abstract0_t*> &, int s, elina_manager_t *man);
    Powerset(std::vector<elina_abstract0_t*>&, std::vector<std::vector<double>>&, int, elina_manager_t *);
    ~Powerset();
    Powerset assign_linexpr_array(
		    elina_dim_t *dims, elina_linexpr0_t** update,
		    unsigned int size,
		    unsigned int output_dim) const;
    Powerset meet_lincons_array(
		    elina_lincons0_array_t* cons) const;
    Powerset join(const Powerset &other) const;

    bool is_bottom() const;
    ssize_t dims() const;
    Powerset &operator=(const Powerset &other);
};

Powerset::Powerset(const Powerset& other) {
  man = other.man;

  disjuncts = std::vector<elina_abstract0_t*>();
  for (elina_abstract0_t* it : other.disjuncts) {
    disjuncts.push_back(elina_abstract0_copy(man, it));
  }

  size = other.size;
  centers = std::vector<std::vector<double>>(other.centers);
}

Powerset::Powerset(std::vector<elina_abstract0_t*>& ds,
    std::vector<std::vector<double>>& cs, int s, elina_manager_t* m) {
  size = s;
  man = m;

  disjuncts = std::vector<elina_abstract0_t*>();
  for (elina_abstract0_t* it : ds) {
    disjuncts.push_back(elina_abstract0_copy(man, it));
  }

  centers = std::vector<std::vector<double>>(cs);
}

Powerset::Powerset(std::vector<elina_abstract0_t*>& ds, int s, elina_manager_t* m) {
  size = s;
  man = m;

  disjuncts = std::vector<elina_abstract0_t*>();
  for (elina_abstract0_t* it : ds) {
    disjuncts.push_back(elina_abstract0_copy(man, it));
  }

  centers = std::vector<std::vector<double>>();
  for (elina_abstract0_t* it : disjuncts) {
    centers.push_back(compute_center(man, it));
  }
}

Powerset::~Powerset() {
  for (elina_abstract0_t* it : disjuncts) {
    elina_abstract0_free(man, it);
  }
}

Powerset Powerset::assign_linexpr_array(elina_dim_t* dims,
    elina_linexpr0_t** update, unsigned int size, unsigned int s) const {
  std::vector<elina_abstract0_t*> ds;
  for (elina_abstract0_t* it : this->disjuncts) {

    elina_abstract0_t* it_dim;
    size_t num_dims = elina_abstract0_dimension(man, it).realdim;
    if (s > num_dims) {
      elina_dimchange_t* dc = elina_dimchange_alloc(0, s - num_dims);
      for (unsigned int i = 0; i < s - num_dims; i++) {
        dc->dim[i] = num_dims;
      }
      it_dim = elina_abstract0_add_dimensions(man, false, it, dc, false);
      elina_dimchange_free(dc);
    } else {
      it_dim = elina_abstract0_copy(man, it);
    }
    elina_abstract0_t* abs = elina_abstract0_assign_linexpr_array(
        this->man, false, it_dim, dims, update, size, NULL);
    if (num_dims > s) {
      elina_dimchange_t* dc = elina_dimchange_alloc(0, num_dims - s);
      for (unsigned int i = 0; i < num_dims - s; i++) {
        dc->dim[i] = s + i;
      }
      abs = elina_abstract0_remove_dimensions(man, true, abs, dc);
      elina_dimchange_free(dc);
    }
    // If the new disjunct is bottom there's no need to add it to the powerset,
    // so we'll drop it for efficiency.
    bool bot = elina_abstract0_is_bottom(this->man, abs);
    if (!bot) {
      ds.push_back(abs);
    } else {
      elina_abstract0_free(this->man, abs);
    }
    elina_abstract0_free(man, it_dim);
  }
  Powerset p(ds, this->size, this->man);
  for (elina_abstract0_t* it : ds) {
    elina_abstract0_free(this->man, it);
  }
  return p;
}

Powerset Powerset::meet_lincons_array(elina_lincons0_array_t* cons) const {
  std::vector<elina_abstract0_t*> ds;
  for (elina_abstract0_t* it : this->disjuncts) {
    elina_abstract0_t* abs = elina_abstract0_meet_lincons_array(
        this->man, false, it, cons);
    bool bot = elina_abstract0_is_bottom(this->man, abs);
    if (!bot) {
      ds.push_back(abs);
    } else {
      elina_abstract0_free(this->man, abs);
    }
  }
  Powerset p(ds, this->size, this->man);
  for (elina_abstract0_t* it : ds) {
    elina_abstract0_free(this->man, it);
  }
  return p;
}

// NOTE: by convention, the this and other should point to the same manager.
// The manager of this is used, so if it is not compatible with the manager of
// other I'm not sure what happens.
Powerset Powerset::join(const Powerset& other) const {
  if (this->disjuncts.size() == 0) {
    return other;
  } else if (other.disjuncts.size() == 0) {
    return *this;
  }
  std::vector<elina_abstract0_t*> ds;
  for (elina_abstract0_t* it : this->disjuncts) {
    ds.push_back(elina_abstract0_copy(man, it));
  }
  for (elina_abstract0_t* it : other.disjuncts) {
    ds.push_back(elina_abstract0_copy(man, it));
  }
  std::vector<std::vector<double>> cs(centers);
  cs.insert(cs.end(), other.centers.begin(), other.centers.end());
  // Now all holds all of the disjuncts from both vectors, so we need to join
  // individual elements until the total number of disjuncts is amll enough
  unsigned int s = std::max(size, other.size);

  while (ds.size() > s) {
    int best_i = 0, best_j = 1;
    double best_dist = norm(cs[0], cs[1]);
    for (unsigned int i = 0; i < ds.size(); i++) {
      for (unsigned int j = i+1; j < ds.size(); j++) {
        double dist = norm(cs[i], cs[j]);
        if (dist < best_dist) {
          best_i = i;
          best_j = j;
          best_dist = dist;
        }
      }
    }
    elina_abstract0_t* n = elina_abstract0_join(man, false, ds[best_i], ds[best_j]);
    elina_abstract0_free(man, ds[best_i]);
    elina_abstract0_free(man, ds[best_j]);
    // j > i so we don't need to worry about messing up indices if we erase
    // j first
    ds.erase(ds.begin() + best_j);
    ds.erase(ds.begin() + best_i);
    cs.erase(cs.begin() + best_j);
    cs.erase(cs.begin() + best_i);
    ds.push_back(n);
    cs.push_back(compute_center(man, n));
  }

  Powerset p(ds, cs, s, man);
  for (elina_abstract0_t* it : ds) {
    elina_abstract0_free(man, it);
  }
  return p;
}

bool Powerset::is_bottom() const {
  bool bottom = true;
  for (elina_abstract0_t* it : this->disjuncts) {
    if (!elina_abstract0_is_bottom(man, it)) {
      bottom = false;
      break;
    }
  }
  // Note that if disjuncts is empty then this powerset is bottom.

  return bottom;
}

Powerset& Powerset::operator=(const Powerset& other) {
  this->size = other.size;
  this->man = other.man;

  // If we're going to change this powerset, we need to free our disjuncts
  // first
  for (elina_abstract0_t* it : this->disjuncts) {
    elina_abstract0_free(this->man, it);
  }

  this->disjuncts = std::vector<elina_abstract0_t*>();
  for (elina_abstract0_t* it : other.disjuncts) {
    this->disjuncts.push_back(elina_abstract0_copy(man, it));
  }

  centers = std::vector<std::vector<double>>(other.centers);
  return *this;
}

typedef struct _bounds_t {
  std::vector<double> lower;
  std::vector<double> upper;
} Bounds;

typedef struct _affineTransform_t {
  std::vector<std::vector<double>> weights;
  std::vector<double> biases;
  int inputSize;
  int outputSize;
} AffineTransform;


std::vector<Powerset> propagatePowerset(const Powerset& p, std::vector<AffineTransform> &t) {
  Powerset z = p;
  std::vector<Powerset> ret;
  ret.push_back(z);
  int num_layers = t.size();
  for (int i = 0; i < num_layers; i++) {
    // Affine transformation
    int in_size = t[i].inputSize;
    int out_size = t[i].outputSize;
    elina_dim_t* dims = (elina_dim_t*) malloc(out_size * sizeof(elina_dim_t));
    elina_linexpr0_t** update = (elina_linexpr0_t**) malloc(out_size *
        sizeof(elina_linexpr0_t*));
    for (int j = 0; j < out_size; j++) {
      dims[j] = j;
      update[j] = elina_linexpr0_alloc(ELINA_LINEXPR_DENSE, in_size);
      for (int k = 0; k < in_size; k++) {
        elina_linexpr0_set_coeff_scalar_double(update[j], k, t[i].weights[j][k]);
      }
      elina_linexpr0_set_cst_scalar_double(update[j], t[i].biases[j]);
    }

    z = z.assign_linexpr_array(dims, update, out_size, out_size);
   
    free(dims);
    for (int j = 0; j < out_size; j++) {
      elina_linexpr0_free(update[j]);
    }
    free(update);
    ret.push_back(z);
    std::cout << "Layer: " << i << " before ReLU" << std::endl;
    std::cout << z.disjuncts.size() << " disjuncts" << std::endl;
    for (elina_abstract0_t* it : z.disjuncts) {
      elina_lincons0_array_t arr = elina_abstract0_to_lincons_array(p.man, it);
      std::cout << "PRINTING BEFORE RELU in Layer: " << i << std::endl;
      elina_lincons0_array_fprint(stdout, &arr, NULL);
      std::cout << "PRINTED BEFORE RELU in Layer: " << i << std::endl;
      elina_lincons0_array_clear(&arr);
    }

    if (i < num_layers - 1) {
      // ReLU
      for (int j = 0; j < out_size; j++) {
        // Match >0 constraint independently along each dimension
        elina_linexpr0_t* lt0_le = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 1);
        elina_linexpr0_t* gt0_le = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 1);
        elina_linexpr0_set_coeff_scalar_double(lt0_le, j, -1.0);
        elina_linexpr0_set_coeff_scalar_double(gt0_le, j, 1.0);

        elina_lincons0_array_t lt0 = elina_lincons0_array_make(1);
        elina_lincons0_array_t gt0 = elina_lincons0_array_make(1);
        lt0.p[0].constyp = ELINA_CONS_SUPEQ;
        lt0.p[0].linexpr0 = lt0_le;
        gt0.p[0].constyp = ELINA_CONS_SUPEQ;
        gt0.p[0].linexpr0 = gt0_le;

        // Meet z with x_j <= 0 and x_j >= 0
        Powerset zlt = z.meet_lincons_array(&lt0);
        // Since we don't change anything for the x_j>=0 case we'll put
        // this in z
        z = z.meet_lincons_array(&gt0);

        // Now we need to assign x_j = 0 in zlt
        elina_linexpr0_t* zero = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 0);
        elina_linexpr0_set_cst_scalar_double(zero, 0.0);
        elina_dim_t dim = j;
        zlt = zlt.assign_linexpr_array(&dim, &zero, 1, out_size);

        // Finally join the outputs
        z = z.join(zlt);

        elina_linexpr0_free(zero);
        elina_lincons0_array_clear(&lt0);
        elina_lincons0_array_clear(&gt0);
      }
    }
    std::cout << "Layer: " << i << " after ReLU" << std::endl;
    std::cout << z.disjuncts.size() << " disjuncts" << std::endl;
    for (elina_abstract0_t* it : z.disjuncts) {
      elina_lincons0_array_t arr = elina_abstract0_to_lincons_array(p.man, it);
      std::cout << "PRINTING AFTER RELU in Layer " << i << std::endl;
      elina_lincons0_array_fprint(stdout, &arr, NULL);
      std::cout << "PRINTED AFTER RELU in Layer " << i << std::endl;
      if (shouldClear) {
	std::cout << "Clearing" << std::endl;
        elina_lincons0_array_clear(&arr);
      }
    }
  }

  return ret;
}

Bounds readProperty(std::string filename) {
  Bounds b;	
  std::vector<double> lower;
  std::vector<double> upper;
  std::ifstream in(filename.c_str());
  std::string line;
  while (getline(in, line)) {
    std::istringstream iss(line);
    double l, u;
    iss.get();
    iss >> l;
    iss.get();
    iss.get();
    iss >> u;
    lower.push_back(l);
    upper.push_back(u);
  }
  b.lower = lower;
  b.upper = upper;
  return b;
}

std::vector<double> parseVector(std::string s) {
  // Remove the [] around the vector
  s.erase(0, 1);
  s.erase(s.end() - 1, s.end());

  std::stringstream ss(s);
  std::string tok;
  std::vector<double> elems;
  // Split s on commas
  while (getline(ss, tok, ',')) {
    if (tok[0] == ' ') {
      tok.erase(0, 1);
    }
    elems.push_back(atof(tok.c_str()));
  }
  return elems;
}

std::vector<std::vector<double>> parseMatrix(std::string s) {
  // Remove [[ and ]] from the input string
  s.erase(0, 2);
  s.erase(s.end() - 2, s.end());

  std::stringstream ss(s);
  std::string tok;
  std::vector<std::vector<double>> rows;
  while (getline(ss, tok, ']')) {
    // Erase ", [" from the beginning of tok
    if (tok[0] == ',') {
      tok.erase(0, 3);
    }
    std::vector<double> b = parseVector("[" + tok + "]");
    rows.push_back(b);
  }
  return rows;
}

std::vector<AffineTransform> readTransforms(std::string filename) {
  int num_layers = 0;
  int input_size = 0;
  int output_size = 0;
  std::vector<AffineTransform> ret;
  std::ifstream file(filename.c_str());
  std::string line;
  int i = 0;
  while (getline(file, line)) {
    if (line != "ReLU") {
      throw std::runtime_error("This is not a fully connected network");
    }
    getline(file, line);
    // Now line holds the weight for this layer
    std::vector<std::vector<double>> m = parseMatrix(line);
    input_size = m[i++].size();
    output_size = m.size();
    getline(file, line);
    std::vector<double> b = parseVector(line);
    AffineTransform t = {m, b, input_size, output_size};
    ret.push_back(t);
  }
  return ret; 
}


void printTransforms(std::vector<AffineTransform> &t) {
  for (AffineTransform it : t) {
    for (std::vector<double> jt : it.weights) {
      for (double kt : jt) {
	std::cout << kt << std::endl;
      }
    }	    
  }
}

void printBounds(Bounds &b) {
  std::cout << "b.lower size: " << b.lower.size() << std::endl;
  for (int i = 0; i < b.lower.size(); i++) {
    std::cout << b.lower[i] << b.upper[i] << std::endl;
  }
}

int main(int argc, char **argv) {
  if (argc != 2) {
	std::cout << "Usage: ./a.out [1|0]" << std::endl;
	exit(1);
  }
  shouldClear = atoi(argv[1]);
  std::cout << "Should clear: " << shouldClear << std::endl;
  std::string prop_file = "./input.txt";
  std::string net_file = "./transformations.txt";
  int inputSize;
  std::vector<AffineTransform> transforms = readTransforms(net_file);
  inputSize = transforms[0].inputSize;
  Bounds b = readProperty(prop_file);
  elina_interval_t** e_itv = (elina_interval_t**) malloc(inputSize * sizeof(elina_interval_t*)); 
  for (int i = 0; i < inputSize; i++) {
    e_itv[i] = elina_interval_alloc();
    elina_interval_set_double(e_itv[i], b.lower[i], b.upper[i]);
  }
  elina_manager_t* man = opt_oct_manager_alloc();
  elina_abstract0_t* init = elina_abstract0_of_box(man, 0, inputSize, e_itv);
  std::vector<elina_abstract0_t*> ds;
  ds.push_back(init);
  Powerset input(ds, 1, man);
  elina_abstract0_free(man, init);
  propagatePowerset(input, transforms);
}
