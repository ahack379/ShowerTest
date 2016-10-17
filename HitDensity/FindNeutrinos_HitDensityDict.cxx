// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME FindNeutrinos_HitDensityDict

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
#include "MakeVtx.h"
#include "VtxDensity.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLVtxDensity_Dictionary();
   static void larlitecLcLVtxDensity_TClassManip(TClass*);
   static void *new_larlitecLcLVtxDensity(void *p = 0);
   static void *newArray_larlitecLcLVtxDensity(Long_t size, void *p);
   static void delete_larlitecLcLVtxDensity(void *p);
   static void deleteArray_larlitecLcLVtxDensity(void *p);
   static void destruct_larlitecLcLVtxDensity(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::VtxDensity*)
   {
      ::larlite::VtxDensity *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::VtxDensity));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::VtxDensity", "VtxDensity.h", 25,
                  typeid(::larlite::VtxDensity), DefineBehavior(ptr, ptr),
                  &larlitecLcLVtxDensity_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::VtxDensity) );
      instance.SetNew(&new_larlitecLcLVtxDensity);
      instance.SetNewArray(&newArray_larlitecLcLVtxDensity);
      instance.SetDelete(&delete_larlitecLcLVtxDensity);
      instance.SetDeleteArray(&deleteArray_larlitecLcLVtxDensity);
      instance.SetDestructor(&destruct_larlitecLcLVtxDensity);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::VtxDensity*)
   {
      return GenerateInitInstanceLocal((::larlite::VtxDensity*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::VtxDensity*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLVtxDensity_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::VtxDensity*)0x0)->GetClass();
      larlitecLcLVtxDensity_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLVtxDensity_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecLcLMakeVtx_Dictionary();
   static void larlitecLcLMakeVtx_TClassManip(TClass*);
   static void *new_larlitecLcLMakeVtx(void *p = 0);
   static void *newArray_larlitecLcLMakeVtx(Long_t size, void *p);
   static void delete_larlitecLcLMakeVtx(void *p);
   static void deleteArray_larlitecLcLMakeVtx(void *p);
   static void destruct_larlitecLcLMakeVtx(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::MakeVtx*)
   {
      ::larlite::MakeVtx *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::MakeVtx));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::MakeVtx", "MakeVtx.h", 25,
                  typeid(::larlite::MakeVtx), DefineBehavior(ptr, ptr),
                  &larlitecLcLMakeVtx_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::MakeVtx) );
      instance.SetNew(&new_larlitecLcLMakeVtx);
      instance.SetNewArray(&newArray_larlitecLcLMakeVtx);
      instance.SetDelete(&delete_larlitecLcLMakeVtx);
      instance.SetDeleteArray(&deleteArray_larlitecLcLMakeVtx);
      instance.SetDestructor(&destruct_larlitecLcLMakeVtx);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::MakeVtx*)
   {
      return GenerateInitInstanceLocal((::larlite::MakeVtx*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::MakeVtx*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLMakeVtx_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::MakeVtx*)0x0)->GetClass();
      larlitecLcLMakeVtx_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLMakeVtx_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLVtxDensity(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::VtxDensity : new ::larlite::VtxDensity;
   }
   static void *newArray_larlitecLcLVtxDensity(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::VtxDensity[nElements] : new ::larlite::VtxDensity[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLVtxDensity(void *p) {
      delete ((::larlite::VtxDensity*)p);
   }
   static void deleteArray_larlitecLcLVtxDensity(void *p) {
      delete [] ((::larlite::VtxDensity*)p);
   }
   static void destruct_larlitecLcLVtxDensity(void *p) {
      typedef ::larlite::VtxDensity current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::VtxDensity

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLMakeVtx(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::MakeVtx : new ::larlite::MakeVtx;
   }
   static void *newArray_larlitecLcLMakeVtx(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::MakeVtx[nElements] : new ::larlite::MakeVtx[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLMakeVtx(void *p) {
      delete ((::larlite::MakeVtx*)p);
   }
   static void deleteArray_larlitecLcLMakeVtx(void *p) {
      delete [] ((::larlite::MakeVtx*)p);
   }
   static void destruct_larlitecLcLMakeVtx(void *p) {
      typedef ::larlite::MakeVtx current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::MakeVtx

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
  void TriggerDictionaryInitialization_libFindNeutrinos_HitDensity_Impl() {
    static const char* headers[] = {
"MakeVtx.h",
"VtxDensity.h",
0
    };
    static const char* includePaths[] = {
"/Users/ah673/WorkArea/Root6LArLite/core",
"/usr/local/include",
"/Applications/BasicSoftware/ROOT/root/include",
"/Users/ah673/WorkArea/Root6LArLite/UserDev/FindNeutrinos/HitDensity/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$VtxDensity.h")))  VtxDensity;}
namespace larlite{class __attribute__((annotate("$clingAutoload$MakeVtx.h")))  MakeVtx;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "MakeVtx.h"
#include "VtxDensity.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::MakeVtx", payloadCode, "@",
"larlite::VtxDensity", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libFindNeutrinos_HitDensity",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libFindNeutrinos_HitDensity_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libFindNeutrinos_HitDensity_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libFindNeutrinos_HitDensity() {
  TriggerDictionaryInitialization_libFindNeutrinos_HitDensity_Impl();
}
