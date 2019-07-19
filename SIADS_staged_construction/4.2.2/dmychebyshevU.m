function y=dmychebyshevU(j,t)


switch j
    case 0
        y = 0;
    case 1
        y = 2;
    case 2
        y = 8*t;
    case 3
        y = 24*t.^2 - 4;
    case 4
        y = 64*t.^3 - 24*t;
    case 5
        y = 160*t.^4 - 96*t.^2 + 6;
    case 6
        y = 384*t.^5 - 320*t.^3 + 48*t;
    case 7
        y = 896*t.^6 - 960*t.^4 + 240*t.^2 - 8;
    case 8
        y = 2048*t.^7 - 2688*t.^5 + 960*t.^3 - 80*t;
    case 9
        y = 4608*t.^8 - 7168*t.^6 + 3360*t.^4 - 480*t.^2 + 10;
end

end