module cubic_b_spline

pub fn cubic_b_spline_interp(points int, knots []f64, values []f64, start_deriv f64, end_deriv f64) [][]f64 {
	secdir := new_second_derivative(knots, values, start_deriv, end_deriv)
	x_values := []f64{len: points + 1, cap: points + 1, init: knots[0] +
		((knots.last() - knots[0]) / points) * index}
	mut y_values := []f64{cap: x_values.len}
	mut i := 0
	for x in x_values {
		for i + 1 < knots.len && knots[i + 1] < x {
			i += 1
			println(i)
		}
		y_values << newint__noinv(x, knots[i], knots[i + 1], values[i], values[i + 1],
			secdir[i], secdir[i + 1])
	}
	return [x_values, y_values]
}

pub fn newint_orig(x f64, a f64, b f64, u f64, v f64, up f64, vp f64) f64 {
	/*
	include error handeling for
	assert( b > a );
    assert( x >= a );
	*/
	ba := b - a
	xa := x - a
	inv_ba := 1.0 / ba
	bx := b - x
	ba2 := ba * ba
	lower := (xa * v) + (bx * u)
	c := (xa * xa - ba2) * xa * vp
	d := (bx * bx - ba2) * bx * up
	return (lower + (0.16666666666666666666) * (c + d)) * inv_ba
}

pub fn newint__noinv(x f64, a f64, b f64, u f64, v f64, up f64, vp f64) f64 {
	/*
	include error handeling for
	assert( b > a );
    assert( x >= a );
	*/
	ba := b - a
	xa := x - a
	bx := b - x
	ba2 := ba * ba
	lower := (xa * v) + (bx * u)
	c := (xa * xa - ba2) * xa * vp
	d := (bx * bx - ba2) * bx * up

	return (lower + (0.16666666666666666666) * (c + d)) / ba
}

pub fn newint__vol(x f64, a f64, b f64, u f64, v f64, up f64, vp f64) f64 {
	/*
	include error handeling for
	assert( b > a );
    assert( x >= a );
	*/
	ba := b - a
	xa := x - a
	inv_ba := 1.0 / ba
	bx := b - x
	ba2 := ba * ba
	lower := (xa * v) + (bx * u)
	c := (xa * xa - ba2) * xa * vp
	d := (bx * bx - ba2) * bx * up

	return (lower + (0.16666666666666666666) * (c + d)) * inv_ba
}

pub fn new_second_derivative(knots []f64, values []f64, start_deriv f64, end_deriv f64) []f64 {
	n := knots.len
	mut c_p := []f64{cap: n}
	mut ypp := []f64{cap: n}
	mut new_x := knots[1]
	mut new_y := values[1]
	mut cj := knots[1] - knots[0]
	mut new_dj := (values[1] - values[0]) / cj

	if start_deriv > 0.99e30 {
		c_p << 0.0
		ypp << 0.0
	} else {
		c_p << 0.5
		ypp << (3.0 * (new_dj - start_deriv) / cj)
	}
	mut i := 1
	for i < (n - 1) {
		old_x := new_x
		old_y := new_y

		aj := cj
		old_dj := new_dj

		new_x = knots[i + 1]
		new_y = values[i + 1]

		cj = new_x - old_x
		new_dj = (new_y - old_y) / cj
		bj := 2.0 * (cj + aj)
		inv_denom := 1.0 / (bj - aj * c_p[i - 1])
		dj := 6.0 * (new_dj - old_dj)

		ypp << ((dj - aj * ypp[i - 1]) * inv_denom)
		c_p << (cj * inv_denom)
		i += 1
	}
	if end_deriv > 0.99e30 {
		c_p << 0.0
		ypp << 0.0
	} else {
		// not really sure why these are here but they are in the origional code
		// old_x := new_x
		// old_y := new_y

		aj := cj
		old_dj := new_dj

		cj = 0.0 // this has the same effect as skipping c_n
		new_dj = end_deriv
		bj := 2.0 * (cj + aj)
		inv_denom := 1.0 / (bj - aj * c_p[i - 1])
		dj := 6.0 * (new_dj - old_dj)

		ypp << ((dj - aj * ypp[i - 1]) * inv_denom)
		c_p << (cj * inv_denom)
	}
	for i > 0 {
		i -= 1
		ypp[i] = ypp[i] - c_p[i] * ypp[i + 1]
	}

	return ypp
}
