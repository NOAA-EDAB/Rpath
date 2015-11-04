#Testing for Rcpp functions

library(Rcpp)

cppFunction('
int SplitSetPred(List stanzas, List state){
  int isp, ist, ia, ieco;
  double Bt, pt, Nt;
  
  //stanza parameters
  const int Nsplit             = as<int>(stanzas["Nsplit"]);     
  const NumericVector Nstanzas = as<NumericVector>(stanzas["Nstanzas"]);
  NumericMatrix NageS          = as<NumericMatrix>(stanzas["NageS"]);
  NumericMatrix WageS          = as<NumericMatrix>(stanzas["WageS"]);
  NumericMatrix WWa            = as<NumericMatrix>(stanzas["WWa"]);
  NumericMatrix Age1           = as<NumericMatrix>(stanzas["Age1"]);
  NumericMatrix Age2           = as<NumericMatrix>(stanzas["Age2"]);
  NumericMatrix EcopathCode    = as<NumericMatrix>(stanzas["EcopathCode"]);

  //state parameters
  NumericVector state_BB = as<NumericVector>(state["BB"]);
  NumericVector state_NN = as<NumericVector>(state["NN"]);

  for (isp = 0; isp < Nsplit; isp++){
    for (ist = 0; ist < Nstanzas[isp]; ist++){
      ieco = EcopathCode(isp, ist);
      //Rcpp::Rcout << "ieco " << ieco << std::endl;
      //Rcpp::Rcout << "Age1 " << Age1(isp, ist) << std::endl;
      //Rcpp::Rcout << "Age2 " << Age2(isp, ist) << std::endl;
      //Rcpp::Rcout << " " << std::endl;
      Bt = 1e-30;
      pt = 1e-30;
      Nt = 1e-30;
      for (ia = Age1(isp, ist); ia <= Age2(isp, ist); ia++){
//Rcpp::Rcout << "isp " << isp << std::endl;        
//Rcpp::Rcout << "ist " << ist << std::endl;
//Rcpp::Rcout << "ia " << ia << std::endl;
        Bt = Bt + NageS(ia, isp) * WageS(ia, isp);
        pt = pt + NageS(ia, isp) * WWa(isp, ia);
        Nt = Nt + NageS(ia, isp);
      }
      state_BB[ieco] = Bt;
      state_NN[ieco] = Nt;
    }
  }
  return(0);
}
')

SplitSetPred(REco.init$stanzas, REco.init$start_state)

#Update numbers, weight, and biomass for multistanza groups
deriv <- deriv_vector(REco.init$params, REco.init$start_state, REco.init$forcing,
                      REco.init$fishing, REco.init$stanzas, 1, 1, 0)

cppFunction('
int update_stanzas(List stanzas, List state, List forcing, List deriv, int yr, int mon){
  int isp, ist, ia, ieco, last, first;
  double Su, Gf, Nt;
  int STEPS_PER_YEAR = 12;

  //stanza parameters
  const int Nsplit                   = as<int>(stanzas["Nsplit"]);     
  const NumericVector Nstanzas       = as<NumericVector>(stanzas["Nstanzas"]);
  const NumericVector vBM            = as<NumericVector>(stanzas["vBM"]);
  const NumericVector Wmat           = as<NumericVector>(stanzas["Wmat"]);
  const NumericVector baseEggsStanza = as<NumericVector>(stanzas["baseEggsStanza"]);
  const NumericVector RscaleSplit    = as<NumericVector>(stanzas["RscaleSplit"]);
  const NumericVector RzeroS         = as<NumericVector>(stanzas["RzeroS"]);
  const NumericVector RecPower       = as<NumericVector>(stanzas["RecPower"]);
  const NumericVector vBGFd          = as<NumericVector>(stanzas["vBGFd"]);
  const NumericVector SpawnEnergy    = as<NumericVector>(stanzas["SpawnEnergy"]);
  const NumericVector SpawnX         = as<NumericVector>(stanzas["SpawnX"]);
  const NumericVector baseSpawnBio   = as<NumericVector>(stanzas["baseSpawnBio"]);
  NumericVector SpawnBio             = as<NumericVector>(stanzas["SpawnBio"]);  
  NumericVector EggsStanza           = as<NumericVector>(stanzas["EggsStanza"]);  
  NumericMatrix NageS                = as<NumericMatrix>(stanzas["NageS"]);
  NumericMatrix WageS                = as<NumericMatrix>(stanzas["WageS"]);
  NumericMatrix SplitAlpha           = as<NumericMatrix>(stanzas["SplitAlpha"]);
  NumericMatrix WWa                  = as<NumericMatrix>(stanzas["WWa"]);
  NumericMatrix Age1                 = as<NumericMatrix>(stanzas["Age1"]);
  NumericMatrix Age2                 = as<NumericMatrix>(stanzas["Age2"]);
  NumericMatrix EcopathCode          = as<NumericMatrix>(stanzas["EcopathCode"]);
  NumericMatrix stanzaPred           = as<NumericMatrix>(stanzas["pred"]);

  //state parameters
  const NumericVector state_BB = as<NumericVector>(state["BB"]);
  
  //forcing parameters
  NumericMatrix force_byrecs   = as<NumericMatrix>(forcing["byrecs"]);

  //derivatives
  const NumericVector LossPropToB = as<NumericVector>(deriv["LossPropToB"]);
  const NumericVector FoodGain    = as<NumericVector>(deriv["FoodGain"]);

  for (isp = 1; isp <= Nsplit; isp++){
    // Update numbers and body weights
    SpawnBio[isp] = 0;
    for(ist = 1; ist <= Nstanzas[isp]; ist++){
      ieco = EcopathCode(isp, ist);
      Su = exp(-LossPropToB[ieco] / STEPS_PER_YEAR / state_BB[ieco]);
      Gf = FoodGain[ieco] / stanzaPred[ieco];
      for(ia = Age1(isp, ist); ia <= Age2(isp, ist); ia++){
        NageS(ia, isp) = NageS(ia, isp) * Su;
        WageS(ia, isp) = vBM[isp] * WageS(ia, isp) + Gf * SplitAlpha(ia, isp);
        if(WageS(ia, isp) > Wmat[isp]){
          SpawnBio[isp] += NageS(ia, isp) * (WageS(ia, isp) - Wmat[isp]);
        };
      }
    }
    EggsStanza[isp] = SpawnBio[isp] * SpawnEnergy[isp] * SpawnX[isp] /
                      (SpawnX[isp] - 1.0 + (SpawnBio[isp] / baseSpawnBio[isp]));
    EggsStanza[isp] *= force_byrecs(yr * STEPS_PER_YEAR + mon, ieco);

    // Need to add monthly recruitment

    // now update n and wt looping backward over age
    last  = Age2(isp, Nstanzas[isp]);
    first = Age1(isp, 1);

    Nt = NageS(last, isp) + NageS(last - 1, isp);    
    if(Nt == 0){Nt = 1e-30;}
    
    WageS(last, isp) = (WageS(last, isp) * NageS(last, isp) + WageS(last - 1, isp) * 
                        NageS(last - 1, isp)) / Nt;
    NageS(last, isp) = Nt;

    for(ia = last - 1; ia > first; ia--){
      NageS(ia, isp) = NageS(ia - 1, isp);
      WageS(ia, isp) = WageS(ia - 1, isp);
    }

    //Apply number of eggs to youngest slot. Includes Walter recruit power
    if(baseEggsStanza[isp] > 0){
      NageS(first, isp) = RscaleSplit[isp] * RzeroS[isp] * pow(double(EggsStanza[isp] / 
                          baseEggsStanza[isp]), double(RecPower[isp]));
    }
    WageS(first, isp) = 0;

    //Uses generalized vonB (exponent is d)
    //Added for stability 4/13/07 (Unlucky Friday)
    for(ia = 0; ia <= last; ia++){
      WWa(ia, isp) = pow(double(WageS(ia, isp)), double(vBGFd[isp]));
    }
  }

return(0);
}
')

update_stanzas(REco.init$stanzas, REco.init$start_state, REco.init$forcing, deriv, 
               1, 1)




