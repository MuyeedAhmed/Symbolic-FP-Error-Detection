def test_matrix_multiplication_errors_on_shapes_mismatch(a, b) -> None:
    t = table_from_pandas(pd.DataFrame({"a": [a], "b": [b]}))
    t.select(c=t.a @ t.b)
    with pytest.raises(ValueError):
        run_all()