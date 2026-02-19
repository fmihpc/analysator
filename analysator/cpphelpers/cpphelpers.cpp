#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <cstdint>
#include <iostream>
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>
#include <sstream>
#include <stdint.h>
#include <unordered_map>
#include <vector>
using namespace std;

// C's % operator is "remainder" operator not modulus like Python's % (not sure
// if there is implementation of this in stnadard library, but its not big so
// here it is)
constexpr int mod(int a,int b) noexcept { return ((a % b) + b) % b; }
constexpr int CHILDS = 8;
constexpr int floordiv(int a,int b) noexcept {
  int q = a / b;
  int r = a % b;
  if ((r != 0) && ((r > 0) != (b > 0))) {
    q--;
  }
  return q;
}

//Convert unordered_map into a Python dictionary
static PyObject *convertToDict(const unordered_map<int, uint64_t> &map) {
  PyObject *dict = PyDict_New();
  for (const auto &it : map) { //structured binding?
    PyObject *key = PyLong_FromLongLong(it.first);
    PyObject *val = PyLong_FromLongLong(it.second);
    PyDict_SetItem(dict, key, val); //error code handling?
    Py_DECREF(key);
    Py_DECREF(val);
  }
  return dict;
}

//could be changed to span?
static void children(int cid, int level, const vector<int64_t> &cid_offsets,
                     const vector<int64_t> &xcells, const vector<int64_t> &ycells,
                     const vector<int64_t> &zcells, vector<int64_t> &out,
                     const vector<vector<int32_t>>& delta
                     ) {
  long cellid = cid - 1 - cid_offsets[level];
  vector<int32_t> cellind(3, -1); //could be reused
  cellind[0] = mod(cellid, (xcells[level])) * 2;
  cellind[1] = mod(floordiv(cellid, xcells[level]), (ycells[level])) * 2;
  cellind[2] = floordiv(cellid, xcells[level] * ycells[level]) * 2;
  //children vector always size 8 (hopefully)
  for (size_t i = 0; i < CHILDS; ++i) {
    out[i] =
        cid_offsets[level + 1] + (cellind[0] + delta[i][0]) +
        xcells[level + 1] * (cellind[1] + delta[i][1]) +
        (cellind[2] + delta[i][2]) * xcells[level + 1] * ycells[level + 1] + 1;
  }

  // return out;
}

//Convert python dictionary to unordered_map
static int convertToUnordMap(PyObject *dict,
                             unordered_map<int, uint64_t> &map) {
  PyObject *key, *value;

  Py_ssize_t pos = 0;

  while (PyDict_Next(dict, &pos, &key, &value)) {

    uint64_t val = PyLong_AsUnsignedLongLong(value);
    long keyval = PyLong_AsLong(key);
    map[keyval] = val;
    if (PyErr_Occurred()) {
      return 1;
    }
  }
  return 0;
}

static PyObject *pyBuildDescriptor(PyObject *self, PyObject *args) {
  unsigned int max_ref_level;
  PyObject *fileindex_for_cellid;
  unsigned int xc, yc, zc;

  stringstream descr; //would vector with push_back be faster?
  // O!|OO (1 required arg (PythonObject) with 2 optional (not sure why we need
  // ! on the first one))
  //  more args O!|O|i for integer
  if (!PyArg_ParseTuple(args, "O|i|i|i|i", &fileindex_for_cellid, &xc, &yc, &zc,
                        &max_ref_level)) {
    return NULL;
  }

  vector<int64_t> xcells(max_ref_level + 1, 0);
  vector<int64_t> ycells(max_ref_level + 1, 0);
  vector<int64_t> zcells(max_ref_level + 1, 0);
  //bitshift very cool
  for (size_t r = 0; r < max_ref_level + 1; r++) {
    xcells[r] = xc << r;
    ycells[r] = yc << r;
    zcells[r] = zc << r;
  }
  unordered_map<int, uint64_t> idxToFileIndex;
  unordered_map<int, uint64_t> fileindex_for_cellid_map;
  if (convertToUnordMap(fileindex_for_cellid, fileindex_for_cellid_map) != 0) {
    cout << "error" << endl;
  }
  int idx = 0;
  vector<vector<int>> subdivided(max_ref_level + 1);
  vector<int64_t> cid_offsets(max_ref_level + 1, 0);
  uint64_t isum = 0;
  for (size_t p = 0; p < max_ref_level; p++) {
    isum = isum + pow(2, 3 * p) * xc * yc * zc;
    cid_offsets[p + 1] = isum;
  }
  for (size_t c = 1; c < xc * yc * zc + 1; c++) {
    if (fileindex_for_cellid_map.find(c) != fileindex_for_cellid_map.end()) {
      // Write
      descr.put('.');
      idxToFileIndex[idx] = fileindex_for_cellid_map[c];
    } else {
      descr.put('R');
      subdivided[0].push_back(c);
    }
    idx = idx + 1;
  }
  static const vector<vector<int32_t>> delta = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
                                   {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}; //array preferable?

  vector<int64_t> childs(8, 0); //array preferabl?
  descr.put('|');
  for (size_t l = 1; l < max_ref_level + 1; l++) {
    auto &vecptr = subdivided[l - 1];
    auto &subd = subdivided[l];
    for (int it : vecptr) {
      children(it, l - 1, cid_offsets, xcells, ycells, zcells, childs, delta);
      for (int64_t child : childs) {
        auto it2 = fileindex_for_cellid_map.find(child);
        if (it2 != fileindex_for_cellid_map.end()) {
          descr.put('.');
          idxToFileIndex[idx] = it2->second;
        } else {
          descr.put('R');
          subd.push_back(child);
        }
        idx = idx + 1;
      }
    }
    if (l < max_ref_level) {
      descr.put('|');
    }
  }
  PyObject *outdict = convertToDict(idxToFileIndex);
  return Py_BuildValue("sO", descr.str().data(), outdict);
}

static PyMethodDef cpphelpers_methods[] = {
    {"buildDescriptor", pyBuildDescriptor, METH_VARARGS, "Build descriptor"}, {NULL} /* Sentinel */
};
static struct PyModuleDef cpphelpers = {
    .m_methods = cpphelpers_methods,
};
PyMODINIT_FUNC PyInit_cpphelpers(void)
// create the module
{
  import_array();
  return PyModuleDef_Init(&cpphelpers);
}
