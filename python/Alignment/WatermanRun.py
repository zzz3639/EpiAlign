import WatermanFun as WF

a=[1,2,3,4,2]
b=[1,2,3,4,3,2,1,2,3,4]
matrix=WF.WatermanAligner_Even(a,b,WF.matchfun_naive,1)
score, matches=WF.TraceBack_Even(a,b,matrix,2)
matches