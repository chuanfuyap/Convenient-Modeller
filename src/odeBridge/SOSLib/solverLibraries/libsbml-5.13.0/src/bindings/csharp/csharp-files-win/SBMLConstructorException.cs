/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 2.0.6
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

namespace libsbml {

 using System;
 using System.Runtime.InteropServices;

/** 
 * @sbmlpackage{core}
 *
@htmlinclude pkg-marker-core.html Exceptions thrown by some libSBML constructors.
 *
 * In some situations, constructors for SBML objects may need to indicate to
 * callers that the creation of the object failed.  The failure may be for
 * different reasons, such as an attempt to use invalid parameters or a
 * system condition such as a memory error.  To communicate this to callers,
 * those classes will throw an SBMLConstructorException.
 *
 * In languages that don't have an exception mechanism (e.g., C), the
 * constructors generally try to return an error code instead of throwing
 * an exception.
 */

public class SBMLConstructorException : System.ArgumentException, IDisposable {
	private HandleRef swigCPtr;
	protected bool swigCMemOwn;
	
	internal SBMLConstructorException(IntPtr cPtr, bool cMemoryOwn)
	{
		swigCMemOwn = cMemoryOwn;
		swigCPtr    = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(SBMLConstructorException obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (SBMLConstructorException obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~SBMLConstructorException() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_SBMLConstructorException(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
    }
  }

  internal SBMLConstructorException(IntPtr cPtr, bool cMemoryOwn, string v) : base(v)
  {
    swigCMemOwn = cMemoryOwn;
    swigCPtr    = new HandleRef(this, cPtr);
  }

  public SBMLConstructorException(string v) : 
   this(libsbmlPINVOKE.new_SBMLConstructorException__SWIG_0(), true, v) 
  {}

  
/** */ /* libsbml-internal */ public
 SBMLConstructorException() : this(libsbmlPINVOKE.new_SBMLConstructorException__SWIG_0(), true) {
  }

  
/** */ /* libsbml-internal */ public
 SBMLConstructorException(string errmsg, string sbmlErrMsg) : this(libsbmlPINVOKE.new_SBMLConstructorException__SWIG_1(errmsg, sbmlErrMsg), true) {
  }

  
/** */ /* libsbml-internal */ public
 SBMLConstructorException(string elementName, SBMLNamespaces xmlns) : this(libsbmlPINVOKE.new_SBMLConstructorException__SWIG_2(elementName, SBMLNamespaces.getCPtr(xmlns)), true) {
  }

  
/**
   * Returns the message associated with this SBML exception.
   *
   * @return the message string.
   */ public
 string getSBMLErrMsg() {
    string ret = libsbmlPINVOKE.SBMLConstructorException_getSBMLErrMsg(swigCPtr);
    return ret;
  }

}

}
