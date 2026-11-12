Exceptions
==========

.. autoclass:: malva_client.exceptions.MalvaAPIError
   :members:
   :show-inheritance:

.. autoclass:: malva_client.exceptions.AuthenticationError
   :members:
   :show-inheritance:

.. autoclass:: malva_client.exceptions.ClientVersionError
   :members:
   :show-inheritance:

.. autoclass:: malva_client.exceptions.SearchError
   :members:
   :show-inheritance:

.. autoclass:: malva_client.exceptions.QuotaExceededError
   :members:
   :show-inheritance:

.. rubric:: Reserved exceptions

``ValidationError``, ``ConfigurationError``, ``StorageError``,
``NetworkError``, and ``TimeoutError`` are kept for compatibility with older
releases. The current client does not raise them.
