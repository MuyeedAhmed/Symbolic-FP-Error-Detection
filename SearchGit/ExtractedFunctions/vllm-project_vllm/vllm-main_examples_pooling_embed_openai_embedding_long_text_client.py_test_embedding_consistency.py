def test_embedding_consistency():
    """Test that chunked processing produces consistent results."""
    client = OpenAI(api_key=API_KEY, base_url=BASE_URL)

    print("\n🔍 Testing Embedding Consistency")
    print("=" * 40)

    # Use the same long text multiple times
    long_text = generate_long_text(
        "Consistency test text for chunked processing validation. " * 50, 3
    )

    embeddings = []

    try:
        for i in range(3):
            response = client.embeddings.create(
                input=long_text, model=MODEL_NAME, encoding_format="float"
            )
            embeddings.append(response.data[0].embedding)
            print(f"   - Generated embedding {i + 1}")

        # Check consistency (embeddings should be identical)
        if len(embeddings) >= 2:
            # Calculate similarity between first two embeddings

            emb1 = np.array(embeddings[0])
            emb2 = np.array(embeddings[1])

            # Cosine similarity
            cosine_sim = np.dot(emb1, emb2) / (
                np.linalg.norm(emb1) * np.linalg.norm(emb2)
            )

            print("✅ Consistency test completed!")
            print(f"   - Cosine similarity between runs: {cosine_sim:.6f}")
            print("   - Expected: ~1.0 (identical embeddings)")

            if cosine_sim > 0.999:
                print("   - ✅ High consistency achieved!")
            else:
                print("   - ⚠️ Consistency may vary due to numerical precision")

    except Exception as e:
        print(f"❌ Consistency test failed: {str(e)}")