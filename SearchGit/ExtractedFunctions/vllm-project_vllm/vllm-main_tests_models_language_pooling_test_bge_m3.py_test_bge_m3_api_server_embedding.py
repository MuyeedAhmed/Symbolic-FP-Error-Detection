async def test_bge_m3_api_server_embedding(server, pooling_task):
    client = server.get_async_client()

    if pooling_task != "embed":
        with pytest.raises(openai.InternalServerError):
            await run_client_embeddings(
                client,
                MODEL_NAME,
                sentences_1,
            )
        return

    embeddings_list_1 = await run_client_embeddings(
        client,
        MODEL_NAME,
        sentences_1,
    )
    embeddings_list_2 = await run_client_embeddings(
        client,
        MODEL_NAME,
        sentences_2,
    )

    embeddings_1 = torch.tensor(embeddings_list_1)
    embeddings_2 = torch.tensor(embeddings_list_2)
    similarity = embeddings_1 @ embeddings_2.T

    # reference values from BAAI/bge-m3 documentation
    reference = torch.tensor(similarity_reference)

    assert torch.allclose(similarity, reference, rtol=0.01)