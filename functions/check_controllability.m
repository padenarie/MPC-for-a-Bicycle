function rnk = check_controllability(sys)

CO = ctrb(sys.A,sys.B); 
rnk = rank(CO)

end