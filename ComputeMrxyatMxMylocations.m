function [MxyNeg1, MxyPos1, MxyNeg2, MxyPos2] = ComputeMrxyatMxMylocations(Mnpx,Mppx,Mnpy,Mppy,Mx,My)

MxyPos1 = sqrt((abs(Mnpx) + Mx).*(abs(Mnpy) + My));
MxyNeg1 = -MxyPos1;

MxyPos2 = sqrt((abs(Mppx) - Mx).*(abs(Mppy) - My));
MxyNeg2 = -MxyPos2;
end

