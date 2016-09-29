// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME ShowerTest_SumChargeDict

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
#include "CalcCharge.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLCalcCharge_Dictionary();
   static void larlitecLcLCalcCharge_TClassManip(TClass*);
   static void *new_larlitecLcLCalcCharge(void *p = 0);
   static void *newArray_larlitecLcLCalcCharge(Long_t size, void *p);
   static void delete_larlitecLcLCalcCharge(void *p);
   static void deleteArray_larlitecLcLCalcCharge(void *p);
   static void destruct_larlitecLcLCalcCharge(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::CalcCharge*)
   {
      ::larlite::CalcCharge *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::CalcCharge));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::CalcCharge", "CalcCharge.h", 25,
                  typeid(::larlite::CalcCharge), DefineBehavior(ptr, ptr),
                  &larlitecLcLCalcCharge_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::CalcCharge) );
      instance.SetNew(&new_larlitecLcLCalcCharge);
      instance.SetNewArray(&newArray_larlitecLcLCalcCharge);
      instance.SetDelete(&delete_larlitecLcLCalcCharge);
      instance.SetDeleteArray(&deleteArray_larlitecLcLCalcCharge);
      instance.SetDestructor(&destruct_larlitecLcLCalcCharge);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::CalcCharge*)
   {
      return GenerateInitInstanceLocal((::larlite::CalcCharge*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::CalcCharge*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLCalcCharge_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::CalcCharge*)0x0)->GetClass();
      larlitecLcLCalcCharge_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLCalcCharge_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLCalcCharge(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::CalcCharge : new ::larlite::CalcCharge;
   }
   static void *newArray_larlitecLcLCalcCharge(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::CalcCharge[nElements] : new ::larlite::CalcCharge[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLCalcCharge(void *p) {
      delete ((::larlite::CalcCharge*)p);
   }
   static void deleteArray_larlitecLcLCalcCharge(void *p) {
      delete [] ((::larlite::CalcCharge*)p);
   }
   static void destruct_larlitecLcLCalcCharge(void *p) {
      typedef ::larlite::CalcCharge current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::CalcCharge

namespace {
  void TriggerDictionaryInitialization_libShowerTest_SumCharge_Impl() {
    static const char* headers[] = {
"CalcCharge.h",
0
    };
    static const char* includePaths[] = {
"/Users/ah673/WorkArea/Root6LArLite/core",
"/Applications/BasicSoftware/ROOT/root/include",
"/Users/ah673/WorkArea/Root6LArLite/UserDev/ShowerTest/SumCharge/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$CalcCharge.h")))  CalcCharge;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "CalcCharge.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::CalcCharge", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libShowerTest_SumCharge",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libShowerTest_SumCharge_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libShowerTest_SumCharge_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libShowerTest_SumCharge() {
  TriggerDictionaryInitialization_libShowerTest_SumCharge_Impl();
}
