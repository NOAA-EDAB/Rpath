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
int update_stanzas(List stanzas, List state, List deriv, int yr, int mon){
  int isp, ist, ia, ieco;
  double Be, Su, Gf;

  //stanza parameters
  const int Nsplit             = as<int>(stanzas["Nsplit"]);     
  const NumericVector Nstanzas = as<NumericVector>(stanzas["Nstanzas"]);
  const NumericVector vBM      = as<NumericVector>(stanzas["vBM"]);
  const NumericVector Wmat     = as<NumericVector>(stanzas["Wmat"]);
  NumericMatrix NageS          = as<NumericMatrix>(stanzas["NageS"]);
  NumericMatrix WageS          = as<NumericMatrix>(stanzas["WageS"]);
  NumericMatrix SplitAlpha     = as<NumericMatrix>(stanzas["SplitAlpha"]);
  NumericMatrix WWa            = as<NumericMatrix>(stanzas["WWa"]);
  NumericMatrix Age1           = as<NumericMatrix>(stanzas["Age1"]);
  NumericMatrix Age2           = as<NumericMatrix>(stanzas["Age2"]);
  NumericMatrix EcopathCode    = as<NumericMatrix>(stanzas["EcopathCode"]);
  NumericMatrix stanzaPred     = as<NumericMatrix>(stanzas["pred"]);

  //state parameters
  const NumericVector state_BB = as<NumericVector>(state["BB"]);

  //derivatives
  const NumericVector LossPropToB = as<NumericVector>(deriv["LossPropToB"]);
  const NumericVector FoodGain    = as<NumericVector>(deriv["FoodGain"]);

  for (isp = 0; isp < Nsplit; isp++){
    // Update numbers and bidy weights
    Be = 0;
    for(ist = 0; ist < Nstanzas[isp]; ist++){
      ieco = EcopathCode(isp, ist);
      Su = exp(-LossPropToB[ieco] / 12 / state_BB[ieco]);
      Gf = FoodGain[ieco] / stanzaPred[ieco];
      for(ia = Age1(isp, ist); ia <= Age2(isp, ist); ia++){
        NageS(ia, isp) = NageS(ia, isp) * Su;
        WageS(ia, isp) = vBM[isp] * WageS(ia, isp) + Gf * SplitAlpha(ia, isp);
        if(WageS(ia, isp) > Wmat[isp]){
          Be = Be + NageS(ia, isp) * (WageS(ia, isp) - Wmat[isp]);
        };
      }

    }
//Rcpp::Rcout << Nstanzas[isp] << std::endl;    
//Rcpp::Rcout << ieco << std::endl;
}
return(0);
}
')

update_stanzas(REco.init$stanzas, REco.init$start_state, deriv, 1, 1)




