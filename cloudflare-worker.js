// ChemCalc API Worker - Cloudflare Workers
// Serves component database, kij parameters, and DIPPR liquid density data
// Hosted at: api.chemcalc.cn

const COMPS = [
  { id:"H2", f:"H\u2082", cn:"\u6c22\u6c14", en:"Hydrogen", Tc:33.19, Pc:1.2964e6, w:-0.219, MW:2.016, cp:[3.249,4.224e-4,8.301e-8,-6.276e-11,0], cat:"perm", cas:"1333-74-0", Tb:20.28, Vc:0.0000642, lel:4, uel:75, ait:500 },
  { id:"N2", f:"N\u2082", cn:"\u6c2e\u6c14", en:"Nitrogen", Tc:126.20, Pc:3.394e6, w:0.0377, MW:28.014, cp:[3.539,-2.610e-4,7.0e-8,-1.567e-11,0], cat:"perm", cas:"7727-37-9", Tb:77.35, Vc:0.0000892 },
  { id:"O2", f:"O\u2082", cn:"\u6c27\u6c14", en:"Oxygen", Tc:154.58, Pc:5.043e6, w:0.0222, MW:31.998, cp:[3.630,-1.794e-3,6.581e-6,-6.014e-9,1.762e-12], cat:"perm", cas:"7782-44-7", Tb:90.19, Vc:0.0000734 },
  { id:"Ar", f:"Ar", cn:"\u6c29\u6c14", en:"Argon", Tc:150.86, Pc:4.898e6, w:0.0, MW:39.948, cp:[2.5,0,0,0,0], cat:"perm", cas:"7440-37-1", Tb:87.3, Vc:0.0000746 },
  { id:"He", f:"He", cn:"\u6c26\u6c14", en:"Helium", Tc:5.19, Pc:0.227e6, w:-0.39, MW:4.003, cp:[2.5,0,0,0,0], cat:"perm", cas:"7440-59-7", Tb:4.22, Vc:0.0000575 },
  { id:"Ne", f:"Ne", cn:"\u6c16\u6c14", en:"Neon", Tc:44.49, Pc:2.679e6, w:-0.029, MW:20.18, cp:[2.5,0,0,0,0], cat:"perm", cas:"7440-01-9", Tb:27.07, Vc:0.0000417 },
  { id:"CO", f:"CO", cn:"\u4e00\u6c27\u5316\u78b3", en:"Carbon Monoxide", Tc:132.92, Pc:3.499e6, w:0.0482, MW:28.01, cp:[3.710,-1.619e-3,3.692e-6,-2.032e-9,2.395e-13], cat:"syn", cas:"630-08-0", Tb:81.66, Vc:0.0000931, lel:12.5, uel:74, ait:609 },
  { id:"CO2", f:"CO\u2082", cn:"\u4e8c\u6c27\u5316\u78b3", en:"Carbon Dioxide", Tc:304.21, Pc:7.383e6, w:0.2236, MW:44.01, cp:[2.401,8.735e-3,-6.607e-6,2.002e-9,0], cat:"syn", cas:"124-38-9", Tb:194.65, Vc:0.000094 },
  { id:"H2O", f:"H\u2082O", cn:"\u6c34", en:"Water", Tc:647.10, Pc:22.064e6, w:0.3449, MW:18.015, cp:[4.070,-1.108e-3,4.152e-6,-2.964e-9,8.07e-13], cat:"syn", cas:"7732-18-5", Tb:373.15, Vc:0.0000559 },
  { id:"H2S", f:"H\u2082S", cn:"\u786b\u5316\u6c22", en:"Hydrogen Sulfide", Tc:373.53, Pc:8.963e6, w:0.0942, MW:34.082, cp:[4.266,-3.468e-3,1.319e-5,-1.097e-8,3.221e-12], cat:"syn", cas:"7783-06-4", Tb:213.6, Vc:0.0000985, lel:4, uel:44, ait:260 },
  { id:"CH4", f:"CH\u2084", cn:"\u7532\u70f7", en:"Methane", Tc:190.56, Pc:4.599e6, w:0.0115, MW:16.043, cp:[1.702,9.081e-3,-2.164e-6,0,0], cat:"hc", cas:"74-82-8", Tb:111.66, Vc:0.0000986, lel:5, uel:15, ait:537 },
  { id:"C2H4", f:"C\u2082H\u2084", cn:"\u4e59\u70ef", en:"Ethylene", Tc:282.35, Pc:5.042e6, w:0.0862, MW:28.054, cp:[1.424,1.436e-2,-5.681e-6,0,0], cat:"hc", cas:"74-85-1", Tb:169.42, Vc:0.000131, lel:2.7, uel:36, ait:450 },
  { id:"C2H6", f:"C\u2082H\u2086", cn:"\u4e59\u70f7", en:"Ethane", Tc:305.32, Pc:4.872e6, w:0.0995, MW:30.07, cp:[1.131,1.925e-2,-5.561e-6,0,0], cat:"hc", cas:"74-84-0", Tb:184.55, Vc:0.0001455, lel:3, uel:12.4, ait:472 },
  { id:"C3H8", f:"C\u2083H\u2088", cn:"\u4e19\u70f7", en:"Propane", Tc:369.83, Pc:4.248e6, w:0.1523, MW:44.097, cp:[1.213,2.849e-2,-8.824e-6,0,0], cat:"hc", cas:"74-98-6", Tb:231.11, Vc:0.0002, lel:2.1, uel:9.5, ait:450 },
  { id:"C3H6", f:"C\u2083H\u2086", cn:"\u4e19\u70ef", en:"Propylene", Tc:365.57, Pc:4.665e6, w:0.1424, MW:42.081, cp:[1.637,2.224e-2,-7.756e-6,0,0], cat:"hc", cas:"115-07-1", Tb:225.46, Vc:0.000185, lel:2, uel:11.1, ait:455 },
  { id:"nC4", f:"n-C4H10", cn:"\u6b63\u4e01\u70f7", en:"n-Butane", Tc:425.12, Pc:3.796e6, w:0.2002, MW:58.123, cp:[1.935,3.669e-2,-1.138e-5,0,0], cat:"hc", cas:"106-97-8", Tb:272.65, Vc:0.000255, lel:1.8, uel:8.4, ait:405 },
  { id:"iC4", f:"i-C4H10", cn:"\u5f02\u4e01\u70f7", en:"i-Butane", Tc:408.80, Pc:3.640e6, w:0.1835, MW:58.123, cp:[1.677,3.756e-2,-1.2e-5,0,0], cat:"hc", cas:"75-28-5", Tb:261.43, Vc:0.0002627, lel:1.8, uel:8.4, ait:460 },
  { id:"SO2", f:"SO\u2082", cn:"\u4e8c\u6c27\u5316\u786b", en:"Sulfur Dioxide", Tc:430.75, Pc:7.884e6, w:0.2451, MW:64.066, cp:[3.267,5.324e-3,-6.843e-7,-2.103e-9,0], cat:"acid", cas:"7446-09-5", Tb:263.13, Vc:0.000122 },
  { id:"NH3", f:"NH\u2083", cn:"\u6c28", en:"Ammonia", Tc:405.40, Pc:11.353e6, w:0.2526, MW:17.031, cp:[3.578,3.02e-3,-1.86e-7,-6.494e-10,0], cat:"acid", cas:"7664-41-7", Tb:239.82, Vc:0.0000725, lel:15, uel:28, ait:651 },
  { id:"COS", f:"COS", cn:"\u7fb0\u57fa\u786b", en:"Carbonyl Sulfide", Tc:378.77, Pc:6.370e6, w:0.0978, MW:60.075, cp:[2.7,8.1e-3,-4.9e-6,1.2e-9,0], cat:"acid", cas:"463-58-1", Tb:222.87, Vc:0.000137, lel:12, uel:29, ait:250 },
  { id:"MeOH", f:"CH\u2083OH", cn:"\u7532\u9187", en:"Methanol", Tc:512.64, Pc:8.097e6, w:0.5656, MW:32.042, cp:[2.211,1.222e-2,-3.450e-6,0,0], cat:"org", cas:"67-56-1", Tb:337.69, Vc:0.000118, lel:6, uel:36, fp:11, ait:464 },
  { id:"EtOH", f:"C\u2082H\u2085OH", cn:"\u4e59\u9187", en:"Ethanol", Tc:513.92, Pc:6.148e6, w:0.6436, MW:46.069, cp:[1.491,2.091e-2,-6.376e-6,0,0], cat:"org", cas:"64-17-5", Tb:351.44, Vc:0.000167, lel:3.3, uel:19, fp:13, ait:363 },
  { id:"Acetone", f:"(CH\u2083)\u2082CO", cn:"\u4e19\u916e", en:"Acetone", Tc:508.20, Pc:4.700e6, w:0.3065, MW:58.080, cp:[1.712,2.617e-2,-8.467e-6,0,0], cat:"org", cas:"67-64-1", Tb:329.2, Vc:0.000209, lel:2.5, uel:12.8, fp:-20, ait:465 },
  { id:"HCHO", f:"CH\u2082O", cn:"\u7532\u919b", en:"Formaldehyde", Tc:408.00, Pc:6.590e6, w:0.2818, MW:30.026, cp:[1.962,1.160e-2,-3.640e-6,0,0], cat:"org", cas:"50-00-0", Tb:254.05, Vc:0.000115, lel:7, uel:73, ait:424 },
  { id:"Benzene", f:"C\u2086H\u2086", cn:"\u82ef", en:"Benzene", Tc:562.05, Pc:4.895e6, w:0.2103, MW:78.114, cp:[-0.206,3.926e-2,-1.330e-5,0,0], cat:"org", cas:"71-43-2", Tb:353.24, Vc:0.000259, lel:1.2, uel:7.8, fp:-11, ait:498 },
  { id:"Toluene", f:"C\u2087H\u2088", cn:"\u7532\u82ef", en:"Toluene", Tc:591.75, Pc:4.108e6, w:0.2640, MW:92.141, cp:[-0.292,4.766e-2,-1.578e-5,0,0], cat:"org", cas:"108-88-3", Tb:383.78, Vc:0.000316, lel:1.1, uel:7.1, fp:4, ait:480 },
  { id:"pXylene", f:"p-C\u2088H\u2081\u2080", cn:"\u5bf9\u4e8c\u7532\u82ef", en:"p-Xylene", Tc:616.20, Pc:3.511e6, w:0.3218, MW:106.167, cp:[-0.286,5.578e-2,-1.846e-5,0,0], cat:"org", cas:"106-42-3", Tb:411.51, Vc:0.000378, lel:1.1, uel:7, fp:27, ait:528 },
  { id:"CyC6", f:"C\u2086H\u2081\u2082", cn:"\u73af\u5df1\u70f7", en:"Cyclohexane", Tc:553.58, Pc:4.075e6, w:0.2094, MW:84.161, cp:[-3.876,6.360e-2,-2.314e-5,0,0], cat:"org", cas:"110-82-7", Tb:353.87, Vc:0.000308, lel:1.3, uel:8, fp:-20, ait:245 },
  { id:"HCl", f:"HCl", cn:"\u6c2f\u5316\u6c22", en:"Hydrogen Chloride", Tc:324.65, Pc:8.310e6, w:0.1319, MW:36.461, cp:[3.512,2.298e-4,-3.560e-7,5.400e-10,0], cat:"acid", cas:"7647-01-0", Tb:188.15, Vc:0.000081 },
  { id:"Cl2", f:"Cl\u2082", cn:"\u6c2f\u6c14", en:"Chlorine", Tc:417.15, Pc:7.711e6, w:0.0688, MW:70.906, cp:[3.266,3.820e-3,-2.230e-6,4.780e-10,0], cat:"acid", cas:"7782-50-5", Tb:239.11, Vc:0.000123 },
  { id:"R22", f:"R22", cn:"R22\u4e8c\u6c1f\u4e00\u6c2f\u7532\u70f7", en:"R22 HCFC-22", Tc:369.30, Pc:4.990e6, w:0.2209, MW:86.468, cp:[2.009,1.810e-2,-7.530e-6,0,0], cat:"ref", cas:"75-45-6", Tb:232.32, Vc:0.000165 },
  { id:"R32", f:"R32", cn:"R32\u4e8c\u6c1f\u7532\u70f7", en:"R32 HFC-32", Tc:351.26, Pc:5.782e6, w:0.2769, MW:52.024, cp:[2.386,1.540e-2,-5.980e-6,0,0], cat:"ref", cas:"75-10-5", Tb:221.5, Vc:0.000121 },
  { id:"R134a", f:"R134a", cn:"R134a\u56db\u6c1f\u4e59\u70f7", en:"R134a HFC-134a", Tc:374.21, Pc:4.059e6, w:0.3268, MW:102.032, cp:[2.014,2.870e-2,-1.220e-5,0,0], cat:"ref", cas:"811-97-2", Tb:247.08, Vc:0.000199 },
  { id:"R125", f:"R125", cn:"R125\u4e94\u6c1f\u4e59\u70f7", en:"R125 HFC-125", Tc:339.17, Pc:3.618e6, w:0.3052, MW:120.022, cp:[2.610,2.920e-2,-1.210e-5,0,0], cat:"ref", cas:"354-33-6", Tb:224.65, Vc:0.00021 },
  { id:"R152a", f:"R152a", cn:"R152a\u4e8c\u6c1f\u4e59\u70f7", en:"R152a HFC-152a", Tc:386.41, Pc:4.517e6, w:0.2752, MW:66.051, cp:[1.803,2.160e-2,-7.870e-6,0,0], cat:"ref", cas:"75-37-6", Tb:249.13, Vc:0.000179 },
  { id:"R1234yf", f:"R1234yf", cn:"R1234yf\u56db\u6c1f\u4e19\u70ef", en:"R1234yf HFO", Tc:367.85, Pc:3.382e6, w:0.276, MW:114.042, cp:[2.500,3.100e-2,-1.350e-5,0,0], cat:"ref", cas:"754-12-1", Tb:243.67, Vc:0.000239 },
];

const KIJ = {"H2-N2":0.0311,"H2-CO":0.0919,"H2-CO2":-0.1622,"H2-CH4":0.0156,"H2-C2H6":0.0089,"H2-C3H8":0.0667,"H2-H2O":0.509,"H2-H2S":0.046,"H2-C2H4":0.0596,"N2-CO":0.0307,"N2-CO2":-0.017,"N2-CH4":0.0311,"N2-C2H6":0.0515,"N2-C3H8":0.0852,"N2-H2O":0.4778,"N2-H2S":0.1672,"N2-O2":-0.0119,"N2-Ar":-0.0026,"N2-NH3":0.2222,"O2-CO2":0.114,"O2-Ar":0.0104,"O2-H2O":0.5,"CO-CO2":-0.054,"CO-CH4":0.03,"CO-H2O":0.49,"CO-H2S":0.0544,"CO2-CH4":0.0919,"CO2-C2H6":0.1322,"CO2-C3H8":0.1241,"CO2-H2O":0.12,"CO2-H2S":0.0974,"CO2-SO2":0.0279,"CO2-COS":0.0394,"CH4-C2H6":-0.0026,"CH4-C2H4":0.0215,"CH4-C3H8":0.014,"CH4-nC4":0.0133,"CH4-iC4":0.0256,"CH4-H2O":0.485,"CH4-H2S":0.08,"C2H6-C3H8":0.0011,"C2H6-H2S":0.0833,"C2H4-C2H6":0.0089,"C3H8-nC4":-0.0078,"C3H8-iC4":-0.0078,"H2S-H2O":0.04,"H2S-COS":0.0259,"NH3-H2O":-0.2589,"nC4-iC4":-0.0004,"CO2-MeOH":0.0504,"CO2-EtOH":0.0896,"CO2-Acetone":0.0478,"CO2-Benzene":0.0774,"CO2-Toluene":0.1033,"H2O-MeOH":-0.0789,"H2O-EtOH":-0.0952,"H2O-Acetone":-0.0608,"MeOH-Benzene":0.09,"MeOH-Toluene":0.085,"Benzene-Toluene":0.0033,"Benzene-pXylene":0.0033,"Toluene-pXylene":0.0011,"Benzene-CyC6":0.0122,"CH4-Benzene":0.0363,"CH4-Toluene":0.0498,"N2-MeOH":0.215,"MeOH-EtOH":0.0,"Acetone-MeOH":0.0,"Acetone-Benzene":0.0,"R32-R125":0.0044,"R134a-R32":0.0,"R134a-R125":0.0,"R22-R134a":0.0,"R22-R125":0.0,"R152a-R134a":0.0};

const LIQDEN = {"N2":{A:3.2091,B:0.2861,C:126.2,D:0.2966},"O2":{A:3.9143,B:0.2861,C:154.58,D:0.2880},"Ar":{A:3.8469,B:0.2881,C:150.86,D:0.2927},"CO":{A:3.2680,B:0.2886,C:132.92,D:0.2780},"CO2":{A:2.768,B:0.2600,C:304.21,D:0.2888},"CH4":{A:2.9214,B:0.2874,C:190.56,D:0.2812},"C2H4":{A:2.5660,B:0.2780,C:282.35,D:0.2864},"C2H6":{A:2.1685,B:0.2723,C:305.32,D:0.2863},"C3H8":{A:1.3757,B:0.2740,C:369.83,D:0.2871},"C3H6":{A:1.4460,B:0.2770,C:365.57,D:0.2920},"nC4":{A:1.0677,B:0.2724,C:425.12,D:0.2863},"iC4":{A:1.0631,B:0.2728,C:408.8,D:0.2790},"H2S":{A:2.7310,B:0.2690,C:373.53,D:0.2710},"SO2":{A:1.9281,B:0.2610,C:430.75,D:0.2780},"NH3":{A:3.5383,B:0.2520,C:405.4,D:0.2860},"COS":{A:1.7904,B:0.2686,C:378.77,D:0.2800},"H2":{A:3.342,B:0.2680,C:33.19,D:0.2765}};

export default {
  async fetch(request) {
    // CORS headers - allow chemcalc.cn and localhost dev
    const origin = request.headers.get('Origin') || '';
    const allowed = ['https://chemcalc.cn', 'https://www.chemcalc.cn', 'http://localhost:5173', 'http://localhost:3000'];
    const corsOrigin = allowed.includes(origin) ? origin : 'https://chemcalc.cn';

    const headers = {
      'Content-Type': 'application/json',
      'Access-Control-Allow-Origin': corsOrigin,
      'Access-Control-Allow-Methods': 'GET, OPTIONS',
      'Access-Control-Allow-Headers': 'Content-Type',
      'Cache-Control': 'public, max-age=3600', // cache 1 hour in browser
    };

    if (request.method === 'OPTIONS') {
      return new Response(null, { status: 204, headers });
    }

    const url = new URL(request.url);

    if (url.pathname === '/api/data') {
      return new Response(JSON.stringify({ COMPS, KIJ, LIQDEN }), { headers });
    }

    return new Response(JSON.stringify({ error: 'Not found' }), { status: 404, headers });
  }
};
