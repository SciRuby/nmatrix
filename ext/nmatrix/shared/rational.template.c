
inline %%TYPE%% BOOL2%%= dtype.id.to_s.upcase%%(bool expr) {
  %%TYPE%% result;
  result.n = expr;
  result.d = 1;
  return result;
}

inline %%TYPE%% %%TYPE_ABBREV%%_bang(%%= dtype.sym == :rational128 ? "int64_t n, int64_t d" : (dtype.sym == :rational64 ? "int32_t n, int32_t d" : "int16_t n, int16_t d")%%)
{
  %%TYPE%% result = {!n, 1};
  return result;
}

inline %%TYPE%% %%TYPE_ABBREV%%_negate(%%= dtype.sym == :rational128 ? "int64_t n, int64_t d" : (dtype.sym == :rational64 ? "int32_t n, int32_t d" : "int16_t n, int16_t d")%%)
{
  %%TYPE%% result = {-n, d};
  return result;
}

inline %%TYPE%% %%TYPE_ABBREV%%_muldiv(int64_t anum, int64_t aden, int64_t bnum, int64_t bden, char k) {
  %%TYPE%% result;
  int64_t t, g1, g2;

  if (k == '/') { // Switch numerator and denominator for division (and move sign)
    if (bnum < 0) {
      anum = -anum;
      bnum = -bnum;
    }
    t = bnum;
    bnum = bden;
    bden = t;
  }

  g1 = nmrb_gcd(anum, bden);
  g2 = nmrb_gcd(aden, bnum);

  result.n = (anum / g1) * (bnum / g2);
  result.d = (aden / g2) * (bden / g1);

  return result;
}

inline %%TYPE%% %%TYPE_ABBREV%%_addsub(int64_t anum, int64_t aden, int64_t bnum, int64_t bden, char k) {
  %%TYPE%% result;

  int64_t ig = nmrb_gcd(aden, bden);
  int64_t a  = anum * (bden / ig);
  int64_t b  = bnum * (aden / ig);
  int64_t c;

  if (k == '+') c=a+b;
  else          c=a-b;

  b        = aden / ig;
  ig       = nmrb_gcd(aden, ig);
  result.n = c / ig;
  a        = bden / ig;
  result.d = a*b;

  return result;
}

inline %%TYPE%% %%TYPE_ABBREV%%_mod(%%= dtype.sym == :rational128 ? "int64_t anum, int64_t aden, int64_t bnum, int64_t bden" : (dtype.sym == :rational64 ? "int32_t anum, int32_t aden, int32_t bnum, int32_t bden" : "int16_t anum, int16_t aden, int16_t bnum, int16_t bden")%%)
{
  // a - (b * int(a/b))
  return %%TYPE_ABBREV%%_addsub(anum, aden, bnum*((int64_t)((anum * bden) / (aden * bnum))), bden, '-');
}
