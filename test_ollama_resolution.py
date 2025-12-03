#!/usr/bin/env python3
"""
Test script for Ollama-based tool resolution
Demonstrates that local AI models can identify EMBOSS tools
"""

import sys
sys.path.insert(0, 'src')

from emboss_wrapper import EMBOSSWrapper


def test_ollama_resolution():
    """Test Ollama's ability to resolve tool names"""
    
    print("=" * 60)
    print("Testing Ollama Tool Resolution")
    print("=" * 60)
    
    # Initialize with Ollama (local mode)
    print("\n🔧 Initializing EMBOSS wrapper with Ollama...")
    emboss = EMBOSSWrapper(ai_mode="local", ai_model="llama3.2:latest")
    
    # Test cases - tools NOT in the explicit tool_map
    test_cases = [
        ("stretcher", "Global alignment tool for longer sequences"),
        ("matcher", "Local alignment showing all matches"),
        ("newcpgseek", "Find CpG islands"),
        ("pepnet", "Protein hydrophobicity plot"),
        ("antigenic", "Predict protein antigenic sites"),
    ]
    
    print("\n📝 Testing tool resolution (these are NOT in tool_map):")
    print("-" * 60)
    
    for query, description in test_cases:
        print(f"\nQuery: '{query}' - {description}")
        resolved = emboss._resolve_tool_with_ai(query)
        
        if resolved == query:
            print(f"  ⚠️  Could not resolve (returned original)")
        else:
            print(f"  ✅ Resolved to: {resolved}")
    
    print("\n" + "=" * 60)
    print("Resolution cache:", emboss.tool_resolution_cache)
    print("=" * 60)


def test_gemini_resolution():
    """Test Gemini's ability to resolve tool names for comparison"""
    
    print("\n\n" + "=" * 60)
    print("Testing Gemini Tool Resolution (for comparison)")
    print("=" * 60)
    
    # Initialize with Gemini (cloud mode)
    print("\n🔧 Initializing EMBOSS wrapper with Gemini...")
    emboss = EMBOSSWrapper(ai_mode="cloud", ai_model="gemini-2.0-flash-exp")
    
    # Same test cases
    test_cases = [
        ("stretcher", "Global alignment tool for longer sequences"),
        ("matcher", "Local alignment showing all matches"),
    ]
    
    print("\n📝 Testing tool resolution:")
    print("-" * 60)
    
    for query, description in test_cases:
        print(f"\nQuery: '{query}' - {description}")
        resolved = emboss._resolve_tool_with_ai(query)
        
        if resolved == query:
            print(f"  ⚠️  Could not resolve (returned original)")
        else:
            print(f"  ✅ Resolved to: {resolved}")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    import os
    
    # Check if Ollama is available
    try:
        import requests
        response = requests.get("http://localhost:11434/api/tags", timeout=2)
        if response.status_code == 200:
            print("✅ Ollama is running")
            test_ollama_resolution()
        else:
            print("❌ Ollama not responding properly")
    except Exception as e:
        print(f"❌ Ollama not available: {e}")
        print("\nTo run this test:")
        print("1. Install Ollama from https://ollama.ai")
        print("2. Run: ollama pull llama3.2")
        print("3. Ollama will auto-start on port 11434")
    
    # Test Gemini if API key available
    if os.getenv('GEMINI_API_KEY'):
        test_gemini_resolution()
    else:
        print("\n💡 Set GEMINI_API_KEY to also test Gemini comparison")
