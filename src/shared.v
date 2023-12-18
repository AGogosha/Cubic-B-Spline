module cubic_b_splinef64

fn newint_orig(x f64, a f64, b f64, u f64, v f64, up f64, vp f64) f64{
	/* 
	include error handeling for 
	assert( b > a );
    assert( x >= a );
	*/
	ba		:= b - a
    xa		:= x - a
    inv_ba	:= 1.0 / ba
    bx		:= b - x
    ba2		:= ba * ba
    lower	:= xa*v + bx*u
    c		:= (xa*xa - ba2)*xa*vp
    d		:= (bx*bx - ba2)*bx*up
	return ( lower + (.16666666666666666666)*( c + d ) ) * inv_ba

}
fn newint__noinv(x f64, a f64, b f64, u f64, v f64, up f64, vp f64) f64{
	/* 
	include error handeling for 
	assert( b > a );
    assert( x >= a );
	*/
    ba		:= b - a
    xa		:= x - a
    bx		:= b - x
    ba2		:= ba * ba
    lower	:= xa*v + bx*u
    c		:= (xa*xa - ba2)*xa*vp
    d		:= (bx*bx - ba2)*bx*up
    
    return ( lower + (.16666666666666666666)*( c + d ) ) / ba
    }

fn newint__vol(x f64, a f64, b f64, u f64, v f64, up f64, vp f64) f64{
	/* 
	include error handeling for 
	assert( b > a );
    assert( x >= a );
	*/
    ba		:= b - a
    xa		:= x - a
	inv_ba	:= 1.0 / ba
    bx		:= b - x
    ba2		:= ba * ba
    lower	:= xa*v + bx*u
    c		:= (xa*xa - ba2)*xa*vp
    d		:= (bx*bx - ba2)*bx*up
    
    return ( lower + (.16666666666666666666)*( c + d ) ) * inv_ba
    }