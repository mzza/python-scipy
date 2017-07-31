def Hamilton_sym(fi):
	''' the del of LEF
		Paremeter:
		----------
		fi: level set fnction
		returns:
		--------
		nabla_sym :a vector,[x,y,z]
	'''
	x=symbols('x',potive=True)
	y=symbols('y',potive=True)
	t=symbols('t',potive=True)
	u=symbols('u',cls=Function)

	u_diff_t=diff(u(x,y,t),t)
	u_diff_x=diff(u(x,y,t),x)
	u_diff_y=diff(u(x,y,t),y)

	fi_diff_t=u_diff_t.subs(u(x,y,t),fi).doit()
	fi_diff_x=u_diff_x.subs(u(x,y,t),fi).doit()
	fi_diff_y=u_diff_y.subs(u(x,y,t),fi).doit()

	nabla_sym=[fi_diff_t,fi_diff_x,fi_diff_y]

	return nabla_sym

def dp_1(s):
	#s=symbols('s',potive=True)
	return sin(2*pi*s)/(2*pi*s)
def dp_2(s):
	#s=symbols('s',potive=True)
	return (s-1)/s

def div_sym(fi_nabla):
	'''
		fi_nabla :vector
	'''
	x=symbols('x',potive=True)
	y=symbols('y',potive=True)
	t=symbols('t',potive=True)
	return diff(fi_nabla[0],t)+diff(fi_nabla[1],x)+diff(fi_nabla[2],y)

def potential_sym_0(fi,mu):
	nabla_sym=Hamilton_sym(fi)
	return mu*div(list(np.array(Hamilton_sym(fi))))

def potential_sym_1(fi,mu):
	nabla_sym=Hamilton_sym(fi)
	len=sqrt(nabla_sym[0]**2+nabla_sym[1]**2+nabla_sym[2]**2)
	return mu*div(list(dp1(len)*np.array(Hamilton_sym(fi))))

def potential_sym_2(fi,mu):
	nabla_sym=Hamilton_sym(fi)
	len=sqrt(nabla_sym[0]**2+nabla_sym[1]**2+nabla_sym[2]**2)
	return mu*div(list(dp2(len)*np.array(Hamilton_sym(fi))))


def potential(fi,mu,t_value,x_value,y_value):
	nabla_sym=Hamilton_sym(fi)
	lam=lambdify([t,x,y],nabla_sym)
	delt=lam(t_value,x_value,y_value)
	len=sqrt(delt[0]**2+delt[1]**2+delt[2]**2)
	if 0<len<1:
		nabla_1=potential_sym_1(fi,mu)
		lam_1=lambdify([t,x,y],nabla_1)
		return lam_1(t_value,x_value,y_value)
	elif len>=1:
		nabla_2=potential_sym_2(fi,mu)
		lam_2=lambdify([t,x,y],nabla_2)
		return lam_2(t_value,x_value,y_value)
	else:
		nabla_0=potential_sym_0(fi,mu)
		lam_2=lambdify([t,x,y],nabla_0)
		return lam_0(t_value,x_value,y_value)
