// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME FindNeutrinos_CalcEfficiencyDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "FilterHitRatio.h"
#include "Sel2CCpi0Eff.h"
#include "SignalEff.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLSignalEff_Dictionary();
   static void larlitecLcLSignalEff_TClassManip(TClass*);
   static void *new_larlitecLcLSignalEff(void *p = 0);
   static void *newArray_larlitecLcLSignalEff(Long_t size, void *p);
   static void delete_larlitecLcLSignalEff(void *p);
   static void deleteArray_larlitecLcLSignalEff(void *p);
   static void destruct_larlitecLcLSignalEff(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::SignalEff*)
   {
      ::larlite::SignalEff *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::SignalEff));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::SignalEff", "SignalEff.h", 25,
                  typeid(::larlite::SignalEff), DefineBehavior(ptr, ptr),
                  &larlitecLcLSignalEff_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::SignalEff) );
      instance.SetNew(&new_larlitecLcLSignalEff);
      instance.SetNewArray(&newArray_larlitecLcLSignalEff);
      instance.SetDelete(&delete_larlitecLcLSignalEff);
      instance.SetDeleteArray(&deleteArray_larlitecLcLSignalEff);
      instance.SetDestructor(&destruct_larlitecLcLSignalEff);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::SignalEff*)
   {
      return GenerateInitInstanceLocal((::larlite::SignalEff*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::SignalEff*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLSignalEff_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::SignalEff*)0x0)->GetClass();
      larlitecLcLSignalEff_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLSignalEff_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecLcLSel2CCpi0Eff_Dictionary();
   static void larlitecLcLSel2CCpi0Eff_TClassManip(TClass*);
   static void *new_larlitecLcLSel2CCpi0Eff(void *p = 0);
   static void *newArray_larlitecLcLSel2CCpi0Eff(Long_t size, void *p);
   static void delete_larlitecLcLSel2CCpi0Eff(void *p);
   static void deleteArray_larlitecLcLSel2CCpi0Eff(void *p);
   static void destruct_larlitecLcLSel2CCpi0Eff(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::Sel2CCpi0Eff*)
   {
      ::larlite::Sel2CCpi0Eff *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::Sel2CCpi0Eff));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::Sel2CCpi0Eff", "Sel2CCpi0Eff.h", 25,
                  typeid(::larlite::Sel2CCpi0Eff), DefineBehavior(ptr, ptr),
                  &larlitecLcLSel2CCpi0Eff_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::Sel2CCpi0Eff) );
      instance.SetNew(&new_larlitecLcLSel2CCpi0Eff);
      instance.SetNewArray(&newArray_larlitecLcLSel2CCpi0Eff);
      instance.SetDelete(&delete_larlitecLcLSel2CCpi0Eff);
      instance.SetDeleteArray(&deleteArray_larlitecLcLSel2CCpi0Eff);
      instance.SetDestructor(&destruct_larlitecLcLSel2CCpi0Eff);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::Sel2CCpi0Eff*)
   {
      return GenerateInitInstanceLocal((::larlite::Sel2CCpi0Eff*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::Sel2CCpi0Eff*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLSel2CCpi0Eff_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::Sel2CCpi0Eff*)0x0)->GetClass();
      larlitecLcLSel2CCpi0Eff_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLSel2CCpi0Eff_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecLcLFilterHitRatio_Dictionary();
   static void larlitecLcLFilterHitRatio_TClassManip(TClass*);
   static void *new_larlitecLcLFilterHitRatio(void *p = 0);
   static void *newArray_larlitecLcLFilterHitRatio(Long_t size, void *p);
   static void delete_larlitecLcLFilterHitRatio(void *p);
   static void deleteArray_larlitecLcLFilterHitRatio(void *p);
   static void destruct_larlitecLcLFilterHitRatio(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::FilterHitRatio*)
   {
      ::larlite::FilterHitRatio *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::FilterHitRatio));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::FilterHitRatio", "FilterHitRatio.h", 25,
                  typeid(::larlite::FilterHitRatio), DefineBehavior(ptr, ptr),
                  &larlitecLcLFilterHitRatio_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::FilterHitRatio) );
      instance.SetNew(&new_larlitecLcLFilterHitRatio);
      instance.SetNewArray(&newArray_larlitecLcLFilterHitRatio);
      instance.SetDelete(&delete_larlitecLcLFilterHitRatio);
      instance.SetDeleteArray(&deleteArray_larlitecLcLFilterHitRatio);
      instance.SetDestructor(&destruct_larlitecLcLFilterHitRatio);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::FilterHitRatio*)
   {
      return GenerateInitInstanceLocal((::larlite::FilterHitRatio*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::FilterHitRatio*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLFilterHitRatio_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::FilterHitRatio*)0x0)->GetClass();
      larlitecLcLFilterHitRatio_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLFilterHitRatio_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLSignalEff(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::SignalEff : new ::larlite::SignalEff;
   }
   static void *newArray_larlitecLcLSignalEff(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::SignalEff[nElements] : new ::larlite::SignalEff[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLSignalEff(void *p) {
      delete ((::larlite::SignalEff*)p);
   }
   static void deleteArray_larlitecLcLSignalEff(void *p) {
      delete [] ((::larlite::SignalEff*)p);
   }
   static void destruct_larlitecLcLSignalEff(void *p) {
      typedef ::larlite::SignalEff current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::SignalEff

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLSel2CCpi0Eff(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::Sel2CCpi0Eff : new ::larlite::Sel2CCpi0Eff;
   }
   static void *newArray_larlitecLcLSel2CCpi0Eff(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::Sel2CCpi0Eff[nElements] : new ::larlite::Sel2CCpi0Eff[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLSel2CCpi0Eff(void *p) {
      delete ((::larlite::Sel2CCpi0Eff*)p);
   }
   static void deleteArray_larlitecLcLSel2CCpi0Eff(void *p) {
      delete [] ((::larlite::Sel2CCpi0Eff*)p);
   }
   static void destruct_larlitecLcLSel2CCpi0Eff(void *p) {
      typedef ::larlite::Sel2CCpi0Eff current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::Sel2CCpi0Eff

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLFilterHitRatio(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::FilterHitRatio : new ::larlite::FilterHitRatio;
   }
   static void *newArray_larlitecLcLFilterHitRatio(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::FilterHitRatio[nElements] : new ::larlite::FilterHitRatio[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLFilterHitRatio(void *p) {
      delete ((::larlite::FilterHitRatio*)p);
   }
   static void deleteArray_larlitecLcLFilterHitRatio(void *p) {
      delete [] ((::larlite::FilterHitRatio*)p);
   }
   static void destruct_larlitecLcLFilterHitRatio(void *p) {
      typedef ::larlite::FilterHitRatio current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::FilterHitRatio

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 457,
                  typeid(vector<int>), DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = 0);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 457,
                  typeid(vector<float>), DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)0x0)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace {
  void TriggerDictionaryInitialization_libFindNeutrinos_CalcEfficiency_Impl() {
    static const char* headers[] = {
"FilterHitRatio.h",
"Sel2CCpi0Eff.h",
"SignalEff.h",
0
    };
    static const char* includePaths[] = {
"/Users/ah673/WorkArea/Root6LArLite/core",
"/usr/local/include",
"/Applications/BasicSoftware/ROOT/root/include",
"/Users/ah673/WorkArea/Root6LArLite/UserDev/FindNeutrinos/CalcEfficiency/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$SignalEff.h")))  SignalEff;}
namespace larlite{class __attribute__((annotate("$clingAutoload$Sel2CCpi0Eff.h")))  Sel2CCpi0Eff;}
namespace larlite{class __attribute__((annotate("$clingAutoload$FilterHitRatio.h")))  FilterHitRatio;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "FilterHitRatio.h"
#include "Sel2CCpi0Eff.h"
#include "SignalEff.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::FilterHitRatio", payloadCode, "@",
"larlite::Sel2CCpi0Eff", payloadCode, "@",
"larlite::SignalEff", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libFindNeutrinos_CalcEfficiency",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libFindNeutrinos_CalcEfficiency_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libFindNeutrinos_CalcEfficiency_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libFindNeutrinos_CalcEfficiency() {
  TriggerDictionaryInitialization_libFindNeutrinos_CalcEfficiency_Impl();
}
