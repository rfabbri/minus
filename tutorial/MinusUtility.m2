needsPackage "MonodromySolver"
debug SLPexpressions -- for newprinttable
-- Utilities to help iterface with Minus
-- To be included ih other files

-- Version of cCode that outputs to file other than stdout ---------------------------------------
cCode (String, List,List) := (outfile, outputs, inputs) -> (
-- template <typename F>
-- inline __attribute__((always_inline)) void 
-- eval<chicago14a, F>::
-- Hxt(const C<F> * __restrict ux, const C<F> * __restrict uparams, C<F> * __restrict uy /*Hxt*/) 
-- header 
    f := openOut outfile;
    h := newPrintTable " = ";
    scan(inputs, g->printName(g,h));
    scan(outputs, g->printName(g,h));
    f << "{ // " << endl;
    scan(h#"#vars",  i-> f << ("  const C<F> &X"|i|" = x["|i|"];") << endl);
    scan(h#"#lines", i-> f << ("  const C<F> "| h#i | ";") << endl);
    scan(#outputs, i-> f << ("  y["|i|"] = "|printName(outputs#i,h)|";") << endl); 
    f << "}";
    close f;
    )

cCode (String, GateMatrix,GateMatrix) := (outfile, M,I) -> cCode(outfile,
flatten entries M, flatten entries I)

-- from original Chicago parser.m2
readStartSys = filename -> (
    l := separate("\n", get filename);
    p0 := value(l#1);
    sols := for i from 3 to #l-2 list value(l#i);
    (transpose matrix p0, sols/(x->transpose matrix x))
    )

-- from original common.m2
-- write starting parameters and solutions to file
writeStartSys = method(Options=>{Filename=>"startSys"})
writeStartSys (Matrix, List) := o -> (M, sols) -> writeStartSys(point M, sols, o)
writeStartSys (Point, List) := o -> (p, sols) -> (
   assert(instance(o.Filename,String));
   f := openOut o.Filename;
   f << "Parameter values: " << endl;
   f << toExternalString p << endl;
   f << "Solutions : " << endl;
   for s in sols do f << toExternalString s << endl;
   close f;
   )
