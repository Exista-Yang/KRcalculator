# Extract selected columns from KnotInfo / LinkInfo xls into JSON (scratch intermediate).
import sys, xlrd, json
src, dst, kind = sys.argv[1], sys.argv[2], sys.argv[3]
wb = xlrd.open_workbook(src, on_demand=True)
sh = wb.sheet_by_index(0)
hdr = [str(c.value) for c in sh.row(0)]
if kind == 'knot':
    keep = ['name','crossing_number','alternating','symmetry_type','pd_notation','braid_notation','braid_index','braid_length',
            'signature','rasmussen_invariant','ozsvath_szabo_tau_invariant','determinant','jones_polynomial','jones_polynomial_vector',
            'homfly_polynomial','homfly_polynomial_vector','khovanov_unreduced_integral_polynomial','khovanov_unreduced_integral_vector',
            'khovanov_reduced_integral_polynomial','khovanov_reduced_integral_vector','khovanov_reduced_rational_vector','dt_notation','three_genus']
else:
    keep = ['name','name_unoriented','knot_atlas','orientation','components','crossing_number','alternating','pd_notation_vector','pd_notation_math',
            'braid_index','braid_notation','braid_notation_old','linking_matrix','signature','nullity','rasmussen_s_invariant','determinant',
            'jones_polynomial','jones_polynomial_vector','homflypt_polynomial','homflypt_polynomial_vector','khovanov_polynomial','khovanov_polynomial_vector',
            'rolfsen_name','dt_code','unoriented']
idx = {h:i for i,h in enumerate(hdr)}
missing = [k for k in keep if k not in idx]
if missing: print('missing columns', missing)
out = []
def cv(c):
    v = c.value
    if c.ctype == xlrd.XL_CELL_NUMBER:
        if float(v).is_integer(): return str(int(v))
        return repr(v)
    return str(v)
for r in range(2, sh.nrows):
    row = sh.row(r)
    d = {k: cv(row[idx[k]]) for k in keep if k in idx}
    out.append(d)
json.dump(out, open(dst,'w'))
print(len(out))
