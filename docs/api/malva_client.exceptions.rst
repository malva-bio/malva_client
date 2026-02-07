Exceptions
==========

All exceptions inherit from :class:`~malva_client.exceptions.MalvaAPIError`,
so you can catch any Malva-related error with a single
``except MalvaAPIError`` clause.

.. autoexception:: malva_client.exceptions.MalvaAPIError
   :members:
   :show-inheritance:

.. autoexception:: malva_client.exceptions.AuthenticationError
   :show-inheritance:

.. autoexception:: malva_client.exceptions.SearchError
   :members:
   :show-inheritance:

.. autoexception:: malva_client.exceptions.QuotaExceededError
   :members:
   :show-inheritance:

.. autoexception:: malva_client.exceptions.ValidationError
   :show-inheritance:

.. autoexception:: malva_client.exceptions.ConfigurationError
   :show-inheritance:

.. autoexception:: malva_client.exceptions.StorageError
   :show-inheritance:

.. autoexception:: malva_client.exceptions.NetworkError
   :show-inheritance:

.. autoexception:: malva_client.exceptions.TimeoutError
   :show-inheritance:
